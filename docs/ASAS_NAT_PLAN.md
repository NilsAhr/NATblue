# ASAS_NAT — Analysis & Implementation Plan

> **Author:** Nils Ahrenhold  
> **Date:** 2026-04-01  
> **Workspace:** NATBlue  
> **Goal:** Implement an independent ASAS (CD&R) optimised for North Atlantic en-route traffic

---

## Table of Contents

1. [Goal 1.1 — Key Components & Functions of the Current ASAS](#11-key-components--functions-of-the-current-asas)
2. [Goal 1.2 — Key Options & Potential Issues for NAT Traffic](#12-key-options--potential-issues-for-nat-traffic)
3. [Goal 1.3 — Solutions & Priority Lists](#13-solutions--priority-lists)
4. [Goal 1.4 — Testing Approach](#14-testing-approach)
5. [Goal 1.5 — Simulation Setup (Scenarios, Logger, Analysis)](#15-simulation-setup)
6. [Goal 1.6 — Testing Strategy & Order](#16-testing-strategy--order)

---

## 1.1 Key Components & Functions of the Current ASAS

### Architecture Overview

```
bluesky/traffic/asas/
├── detection.py      ← ConflictDetection base class (replaceable Entity)
├── statebased.py     ← StateBased CD: matrix-based all-pairs detection (CURRENT, returns 10 values incl. dalt)
├── statebased_old.py ← StateBased CD: original version (returns 9 values, no dalt)
├── resolution.py     ← ConflictResolution base class (replaceable Entity) + resumenav()
├── mvp.py            ← MVP: Modified Voltage Potential resolution (default CR)
└── src_cpp/          ← C++ accelerated detection (returns 9 values — INCOMPATIBLE with current update())
```

```
Update loop (traffic.py):
  1. Atmosphere → 2. ADSB → 3. Autopilot (ap.update)
  4. IF asastimer ready:
       cd.update(traf, traf)       ← Conflict Detection
       cr.update(cd, traf, traf)   ← Conflict Resolution
  5. aporasas.update()             ← Select AP or ASAS commands
  6. Performance limits → 7. Kinematics
```

### ConflictDetection (detection.py)

| Component | Description |
|---|---|
| **Protected zone** | `rpz` (horizontal, NM→m), `hpz` (vertical, ft→m) — per-aircraft arrays |
| **Lookahead** | `dtlookahead` (default 300s = 5 min) — per-aircraft |
| **Outputs** | `confpairs`, `lospairs`, `inconf`, `tcpamax`, `qdr`, `dist`, `dcpa`, `tcpa`, `tLOS`, `dalt` |
| **Unique tracking** | `confpairs_unique`, `lospairs_unique` (frozensets) |
| **Cumulative** | `confpairs_all`, `lospairs_all` (since t=0) |
| **Stack commands** | `CDMETHOD`/`ASAS`, `ZONER`/`PZR`, `ZONEDH`/`PZDH`, `DTLOOK`, `DTNOLOOK` |

**NATBlue modification:** `detection.py` unpacks **10 return values** (added `self.dalt`) in `update()`. The `statebased.py` computes minimum vertical distance during conflict window. The C++ version only returns 9 → **will crash if C++ is selected**.

### StateBased Detection (statebased.py)

Algorithm (all-pairs O(n²) matrix):
1. **Horizontal**: QDR/dist matrix → relative velocity → tCPA → dCPA² → compare to RPZ² → tin_hor/tout_hor
2. **Vertical**: dalt matrix → relative VS → crossing times through ±HPZ → tin_ver/tout_ver
3. **Combined**: tin_conf = max(tin_hor, tin_ver), tout_conf = min(tout_hor, tout_ver) → boolean conflict
4. **dalt_min (NATBlue addition)**: Minimum vertical distance during [tin_conf, tout_conf] — checks entry, exit, and zero-crossing

### ConflictResolution (resolution.py)

| Component | Description |
|---|---|
| **Resolution factors** | `resofach` (1.01 default), `resofacv` (1.01 default) — multiplied on PZ for resolution margin |
| **Per-aircraft** | `active` (bool), `trk`, `tas`, `alt`, `vs` — ASAS-commanded state |
| **Properties** | `hdgactive`, `vsactive`, `altactive`, `tasactive` — all default to `active` |
| **resumenav()** | Recovery logic: past CPA + not in hor-LOS + not bouncing → deactivate ASAS → `route.direct()` to next waypoint |
| **Bouncing check** | `anglediff < 30°` AND `hdist < rpz * resofach` |
| **Stack commands** | `RESO`, `PRIORULES`, `NORESO`, `RESOOFF`, `RFACH`, `RFACV`, `RSZONER`, `RSZONEDH` |

### MVP Resolution (mvp.py)

| Component | Description |
|---|---|
| **Dimension switches** | `swresohoriz` (default True), `swresospd`, `swresohdg`, `swresovert` |
| **Stack commands** | `RMETHH` (BOTH/SPD/HDG/NONE), `RMETHV` (ON/OFF) |
| **Priority codes** | FF1 (no prio), FF2 (cruising prio), FF3 (climbing prio), LAY1/LAY2 (layered) |
| **resolve()** | Loops conflict pairs → calls `MVP()` per pair → accumulates dv vector → applies prio → converts to polar → caps to perf limits |
| **MVP()** | Horizontal: CPA vector, intrusion depth, erratum correction for grazing → dv1, dv2. Vertical: intrusion depth, solve time → dv3 |
| **Speed capping** | `perf.vmin/vmax` for GS, `perf.vsmin/vsmax` for VS |
| **Alt computation** | `vscapped * timesolveV + alt` → compared to `selalt` to avoid overshooting |

### APorASAS (aporasas.py)

Final selector per channel: `trk = where(cr.hdgactive, cr.trk, ap.trk)` etc. Converts ASAS GS→TAS with wind correction. Forces VS positive (sign handled via delalt elsewhere).

---

## 1.2 Key Options & Potential Issues for NAT Traffic

### NAT Traffic Characteristics

| Characteristic | Implication for ASAS |
|---|---|
| **Bidirectional flow** (W→E, E→W) | Head-on (160–180°) and same-direction (0–20°) conflicts dominate |
| **Shallow angles (0–60°)** | Very long conflicts (>10 min), slow relative velocity, CPA far in future |
| **Wide angles (120–180°)** | Fast closing speed, short reaction time, head-on geometry singularities |
| **Parallel tracks** | Aircraft flying on adjacent organized tracks ≈ 30–60 NM apart, lateral offsets small |
| **Continuous climb** | Vertical PZ intrusions during climbs through FL300–FL410 |
| **Near MMO** | Speed buffer nearly zero → MVP cannot use speed-decrease resolutions |
| **Near max altitude** | VS margin nearly zero → vertical resolution capped at ceiling |
| **Long conflicts** | Bouncing prevention must handle 10+ min resolution sequences |
| **14 NM separation (target)** | 2.8× larger PZ than continental (5 NM) → much larger detection/resolution volumes |
| **En-route only** | No terminal manoeuvring, simple straight-line segments between waypoints |

### Issues Identified (Your Analysis)

#### ISSUE-1: Bouncing / Parallel Drift (HIGH RISK for NAT)
- **Symptom**: After resolving, two near-parallel aircraft bounce back into conflict. The bouncing check (`anglediff < 30° AND hdist < rpz*resofach`) keeps ASAS active, but the aircraft drift laterally away from their planned trajectory — potentially into adjacent tracks.
- **Root cause**: MVP provides a dv that pushes aircraft apart, but for shallow-angle conflicts on parallel tracks, this creates oscillations or permanent lateral offset.
- **NAT amplifier**: Organized Track System has parallel tracks — lateral drift can cause secondary conflicts with adjacent traffic.
- **Professor's concern**: "Central sequencing" — aircraft on parallel tracks all resolving simultaneously could create cascading conflicts.

#### ISSUE-2: Speed Resolution Near MMO (HIGH RISK)
- **Symptom**: Aircraft flying close to MMO cannot slow down; `perf.vmax` cap clips the resolution speed.
- **Root cause**: MVP caps `newgs` to `[perf.vmin, perf.vmax]`. For NAT en-route, most aircraft cruise at M0.82–0.85 while MMO ≈ M0.87. The speed margin is tiny.
- **Nino's finding**: GS/TAS conversions create singularities at certain angles (common in NAT headings). He removed unnecessary conversions and used `perf.limits()` for TAS directly.

#### ISSUE-3: Vertical Resolution Near Max Altitude (HIGH RISK)
- **Symptom**: Aircraft at FL390–FL410 cannot climb further; `perf.vsmax` approaches zero near service ceiling.
- **Root cause**: MVP computes `dv3` assuming vertical freedom. Near ceiling, the capped `vscapped` may be zero, making the altitude resolution `asasalttemp = 0 * timesolveV + alt = alt` (no change).
- **Nino's finding**: Autopilot waypoint altitude logic can trigger climbs into vertical PZ of an aircraft in the horizontal PZ.

#### ISSUE-4: Vertical PZ Intrusion During Autopilot Manoeuvres (MEDIUM)
- **Symptom**: Two aircraft within horizontal PZ, one starts a climb at a waypoint → enters the other's vertical PZ.
- **Nino's solution**: VS-stop: predict if conflict pair within horizontal PZ will enter vertical PZ next timestep; halt climb/descent.
- **Not yet implemented in NATBlue.**

#### ISSUE-5: Persistent LoS at Spawning (LOW for well-designed scenarios)
- **Symptom**: Aircraft spawn inside each other's PZ → immediate LoS that MVP struggles with.
- **Nino's solution**: Time-shift initial deconfliction (120s shifts until no interference).
- **Mitigation**: Can be handled in scenario generation (pre-deconfliction) rather than in ASAS.

#### ISSUE-6: Simplistic Recovery Phase Logic (MEDIUM)
- **Symptom**: Current recovery uses only "past CPA" check. For long NAT conflicts, CPA-based recovery may be premature.
- **Nino's solution**: 2-criteria recovery using intruder state at conflict start AND current state.
- **Not yet implemented in NATBlue.**

#### ISSUE-7: RPZ Inconsistency (BUG)
- **Current state**: `detection.py` hardcodes `asas_pzr=1.0`, while `settings.cfg` has `asas_pzr=5.0`. The code default (1.0) overrides the cfg because `set_variable_defaults` only applies if the variable is not already set in code.
- **Fix needed**: Change `detection.py` default to 5.0, then use cfg to set NAT value (14.0).

#### ISSUE-8: Intent Information Not Used
- **Current state**: Detection is state-based only (current position, velocity). No use of FMS route, waypoint, or planned altitude data.
- **Opportunity**: For NAT with known routes, intent-based detection could predict conflicts much earlier and more accurately. Flight plan waypoints are available in `traf.ap.route[idx]`.

---

## 1.3 Solutions & Priority Lists

### Priority 1: Critical Foundation (Must Have)

| # | Solution | Files to Change | Effort |
|---|---|---|---|
| **P1.1** | **Fix RPZ default** — change `detection.py` L10 from `asas_pzr=1.0` to `asas_pzr=5.0` | `detection.py` | 1 line |
| **P1.2** | **Create ASAS_NAT folder** — new MVP_NAT class as a plugin in `bluesky/plugins/asas_nat/` that extends MVP with NAT-specific control logic | New file `mvp_nat.py` | Medium |
| **P1.3** | **Vertical-first resolution logic** — since vertical is most promising for NAT, implement a resolution selector that prefers vertical when feasible, falls back to horizontal | `mvp_nat.py` | Medium |
| **P1.4** | **NAT parameter settings** — configurable via scenario commands: PZR=14NM, HPZ=1000ft, DTLOOK=600s+ | Scenario `.scn` files | Small |

### Priority 2: Core NAT Issues (Should Have)

| # | Solution | Files to Change | Effort |
|---|---|---|---|
| **P2.1** | **VS-stop (ISSUE-4)** — predict vertical PZ entry for aircraft within horizontal PZ; halt climb/descent | `mvp_nat.py` | Medium |
| **P2.2** | **Anti-bouncing for parallel tracks** — enhanced bouncing logic that considers track geometry; limit lateral offset from planned route | `mvp_nat.py` or `resolution.py` override | Medium |
| **P2.3** | **2-criteria recovery (ISSUE-6)** — override `resumenav()` to use both initial and current intruder state | `mvp_nat.py` | Medium |
| **P2.4** | **MMO-aware speed resolution** — before computing speed dv, check available speed margin to MMO; if margin < threshold, suppress speed component and shift to heading or vertical | `mvp_nat.py` | Medium |

### Priority 3: Advanced Enhancements (Nice to Have)

| # | Solution | Files to Change | Effort |
|---|---|---|---|
| **P3.1** | **Intent-based lookahead** — use `traf.ap.route` waypoints to predict future conflicts beyond state-based detection | New detection method or extension | Large |
| **P3.2** | **Altitude-ceiling-aware vertical resolution** — check `perf.hmaxact` before issuing vertical resolution; if near ceiling, fallback to horizontal | `mvp_nat.py` | Medium |
| **P3.3** | **GS/TAS singularity fix (Nino)** — use `perf.limits()` TAS directly, avoid repeated GS↔TAS conversions at NAT headings | `mvp_nat.py` | Small |
| **P3.4** | **Initial deconfliction** — scenario-side time shifting for spawn-in-PZ prevention | Scenario generation script | Medium |

### How to Use Intent Information

Intent information available in BlueSky:
- `traf.ap.route[idx].wplat / wplon / wpalt / wpspd` — full waypoint list
- `traf.ap.route[idx].wpname` — waypoint names
- `traf.ap.route[idx].iactwp` — active waypoint index
- `traf.ap.dest` — destination airport
- `traf.selalt` — selected altitude (FMS target)

**Potential uses:**
1. **Planned trajectory conflict detection** — project aircraft along their route (not just current heading) to detect conflicts 20–30 min ahead. This is especially valuable for crossing traffic where current heading gives poor prediction.
2. **Altitude intent** — if an aircraft is climbing to FL390 (`selalt`), predict when it will be at each FL and check for vertical conflicts with cruising traffic at those levels.
3. **Waypoint-based turn prediction** — predict heading changes at upcoming waypoints to assess if they create new conflicts.
4. **Priority assignment** — aircraft with longer remaining route may get lower priority (more time to recover). Aircraft on continuous climb could get priority (FF3-style).

**Implementation approach:**
- Phase 1: Use `selalt` in vertical resolution logic (check if FMS altitude target resolves the conflict naturally — already partially done in MVP's `signdvs/signalt` logic).
- Phase 2: Build a simple trajectory predictor using waypoint list for lookahead beyond current state.

---

## 1.4 Testing Approach

### What We Need to Verify

| Test Category | Question | Metrics |
|---|---|---|
| **Detection correctness** | Does CD detect all conflicts at 5NM / 14NM? | # conflicts detected, # missed, false positives |
| **Resolution effectiveness** | Does CR prevent LoS? | # LoS events, LoS duration, minimum separation achieved |
| **Resolution quality** | How much does CR disturb the trajectory? | Lateral deviation from plan, altitude deviation, fuel penalty |
| **Recovery** | Do aircraft return to their planned route? | Time to recover, distance from planned waypoint at recovery |
| **Stability** | No cascading / bouncing? | # secondary conflicts, max simultaneous active resolutions |
| **Performance limits** | Does CR respect MMO and ceiling? | Speed/altitude exceedances, failed resolutions |
| **Vertical resolution** | Does vertical-first strategy work? | % of conflicts resolved vertically vs horizontally, LoS rate |

### Test Geometries for NAT

| Geometry | Description | Expected Difficulty |
|---|---|---|
| **G1: Same-direction, same FL** | Two aircraft on same track, different speeds | Easy (speed only) |
| **G2: Same-direction, climbing** | One climbing through other's FL | Medium (VS-stop needed) |
| **G3: Parallel tracks, same FL** | Adjacent tracks, lateral conflict | Hard (bouncing risk) |
| **G4: Parallel tracks, one climbing** | Cross between tracks at different FLs | Hard (combined) |
| **G5: Head-on, same FL** | Opposing traffic on same track | Medium (but fast) |
| **G6: Head-on, offset** | Opposing traffic on adjacent tracks | Medium-Hard |
| **G7: Crossing 30–60°** | Shallow crossing angle | Very Hard (long duration) |
| **G8: Crossing at MMO** | Conflict with both aircraft near MMO | Hard (no speed margin) |
| **G9: Near ceiling** | Conflict at FL410 with limited climb | Hard (no vertical margin) |
| **G10: Multi-aircraft cluster** | 3–4 aircraft in proximity | Very Hard (cascading) |

---

## 1.5 Simulation Setup

### A. Logger Enhancements for ASAS Testing

The current `loggercompact.py` logs good baseline data but needs **ASAS-specific columns** to diagnose resolution behaviour. 

#### Additional Columns Needed in Flight-State Log (flst)

| Column | Source | Purpose |
|---|---|---|
| `asas_active` | `traf.cr.active[idx]` | Is ASAS controlling this aircraft right now? |
| `asas_trk` | `traf.cr.trk[idx]` | ASAS-commanded track (vs. AP track) |
| `asas_tas` | `traf.cr.tas[idx]` | ASAS-commanded speed |
| `asas_vs` | `traf.cr.vs[idx]` | ASAS-commanded vertical speed |
| `asas_alt` | `traf.cr.alt[idx]` | ASAS-commanded altitude |
| `ap_trk` | `traf.ap.trk[idx]` | Autopilot desired track (planned) |
| `ap_alt` | `traf.ap.alt[idx]` | Autopilot desired altitude |
| `selalt` | `traf.selalt[idx]` | FMS selected altitude |
| `mach` | `traf.M[idx]` or compute | Current Mach number |
| `gs` | `traf.gs[idx]` | Ground speed (for GS/TAS analysis) |
| `hdg` | `traf.hdg[idx]` | Actual heading (wind effect) |
| `trk` | `traf.trk[idx]` | Actual track |
| `n_conflicts` | count from `cd.confpairs` | Number of simultaneous conflicts this aircraft is in |
| `lat_deviation` | compute from route | Lateral deviation from planned route [NM] |

#### Additional Columns Needed in Conflict Log (conflog)

| Column | Source | Purpose |
|---|---|---|
| `reso_dimension` | `'V'`/`'H'`/`'HV'` from mvp_nat | Which dimension was the resolution in? |
| `asas_active_ac1` | `traf.cr.active[i1]` | Was ASAS active for ac1? |
| `asas_active_ac2` | `traf.cr.active[i2]` | Was ASAS active for ac2? |
| `mach_ac1` | `traf.M[i1]` | Mach of ac1 (diagnose MMO issues) |
| `mach_ac2` | `traf.M[i2]` | Mach of ac2 |
| `alt_max_ac1` | `traf.perf.hmaxact[i1]` | Max altitude ac1 (diagnose ceiling issues) |
| `selalt_ac1` | `traf.selalt[i1]` | FMS alt target ac1 |
| `selalt_ac2` | `traf.selalt[i2]` | FMS alt target ac2 |
| `conflict_angle` | computed from headings | Relative angle (diagnose geometry) |
| `los_severity` | `min(dist)/rpz` | How deep was the intrusion (0=touch, 1=centre) |
| `resolution_success` | bool | Did this conflict end without LoS? |

### B. Scenario Files

#### B1: Unit Test Scenarios (2 aircraft each)

Create dedicated `.scn` files for each geometry:

```
scenario/ASAS_NAT/
├── test_G1_same_dir_same_FL.scn      ← Same direction, same FL
├── test_G2_same_dir_climbing.scn      ← One climbing through other's FL
├── test_G3_parallel_same_FL.scn       ← Adjacent parallel tracks
├── test_G4_parallel_climbing.scn      ← Parallel, one climbing
├── test_G5_head_on_same_FL.scn        ← Opposing traffic
├── test_G6_head_on_offset.scn         ← Opposing, laterally offset
├── test_G7_crossing_shallow.scn       ← 30-60° crossing
├── test_G8_crossing_MMO.scn           ← Conflict near MMO
├── test_G9_near_ceiling.scn           ← Conflict at FL410
├── test_G10_multi_aircraft.scn        ← 4-aircraft cluster
└── test_ASAS_config.scn               ← Common ASAS settings (included by all)
```

**Common config (test_ASAS_config.scn):**
```
00:00:00.00>CDMETHOD STATEBASED
00:00:00.00>RESO MVP
00:00:00.00>ZONER 5
00:00:00.00>ZONEDH 1000
00:00:00.00>DTLOOK 300
00:00:00.00>RMETHV ON
00:00:00.00>STARTLOG
```

**Example — G5 head-on (test_G5_head_on_same_FL.scn):**
```
# Two aircraft head-on at FL370 mid-Atlantic
00:00:00.00>CALL test_ASAS_config.scn
00:00:00.00>CRE AC001 B77W 53.0 -30.0 090 FL370 480
00:00:00.00>CRE AC002 A333 53.0 -28.0 270 FL370 480
00:00:00.00>DEST AC001 EGLL
00:00:00.00>DEST AC002 KJFK
# Run for 30 minutes
00:00:00.00>FF
00:30:00.00>HOLD
```

#### B2: Stress Test Scenarios (10–20 aircraft)

- Multiple parallel eastbound + westbound streams
- Mix of cruising and climbing aircraft
- Spawn at entry points of a simplified NAT track structure

#### B3: NAT-Border Spawn Scenario (for full-day testing)

**Approach:** Transform the existing `optimised.scn` (1000 flights) into a NAT-border scenario:

1. Parse `optimised.scn` to extract each flight's route
2. Find the waypoint closest to a defined NAT boundary (e.g., 10°W–60°W, 45°N–63°N)
3. For each flight, create a new CRE command at the NAT boundary waypoint with:
   - The flight's speed, altitude, and heading at that point
   - The remaining route from the NAT boundary to destination
4. This dramatically reduces simulation time (skip continental routing)

**Script needed:** `utils/create_nat_border_scenario.py`

### C. Analysis Pipeline

Post-processing should compute from the logs:

| Metric | Formula | Acceptable Value |
|---|---|---|
| **LoS rate** | # LoS events / # conflicts | < 5% (target: 0%) |
| **LoS duration** | Average duration of LoS events [s] | < 30s |
| **LoS severity** | Average `min_dist / rpz` during LoS | > 0.8 (barely touching) |
| **Conflict duration** | Average `toutconf - tinconf` [s] | Informational |
| **Resolution success rate** | # conflicts without LoS / # total conflicts | > 95% |
| **Max lateral deviation** | Max `lat_deviation` during resolution [NM] | < 5 NM (1 track width) |
| **Fuel penalty** | (total_fuel_with_ASAS - total_fuel_baseline) / total_fuel_baseline | < 2% |
| **Secondary conflict rate** | # new conflicts created during resolution / # original conflicts | < 10% |
| **Bouncing events** | # conflicts with duration > 600s that oscillate | 0 target |
| **MMO violations** | # timesteps where M > MMO | 0 |
| **Altitude ceiling violations** | # timesteps where alt > max_alt | 0 |

---

## 1.6 Testing Strategy & Order

### Phase 0: Baseline — Detection Only (No Resolution)

**Goal:** Verify detection works correctly at 5 NM, then at 14 NM. Establish baseline conflict counts.

| Step | Action | Scenario | Duration |
|---|---|---|---|
| 0.1 | Fix RPZ default in `detection.py` (1.0 → 5.0) | — | 5 min |
| 0.2 | Run G5 (head-on) with `RESO OFF` | `test_G5_head_on_same_FL.scn` | 30 min sim |
| 0.3 | Verify conflict detected at correct range, correct tCPA | Manual log check | — |
| 0.4 | Run G1–G10 all with `RESO OFF` | All unit test scenarios | 1 hour |
| 0.5 | Run full-day baseline with `RESO OFF`, PZR=5 then PZR=14 | `optimised.scn` or NAT-border | 2 hours |

**Exit criteria:** All geometries produce expected conflict detections. No crashes, no missed conflicts.

### Phase 1: Basic MVP Vertical Resolution (Start Here)

**Goal:** Test stock MVP with `RMETHV ON` on simple geometries. Identify failure modes.

| Step | Action | Scenario | What to Look For |
|---|---|---|---|
| 1.1 | MVP vertical only on G2 (climbing through) | `test_G2` | Does VS-stop prevent intrusion? |
| 1.2 | MVP vertical only on G5 (head-on same FL) | `test_G5` | Vertical resolution at same FL (should step up/down) |
| 1.3 | MVP vertical only on G9 (near ceiling) | `test_G9` | Does it fail when aircraft can't climb? |
| 1.4 | MVP vertical only on G8 (near MMO) | `test_G8` | Vertical should work regardless of speed |
| 1.5 | MVP vertical only on G3 (parallel tracks) | `test_G3` | Does it oscillate? Does it resolve? |
| 1.6 | MVP vertical on all G1–G10 | All unit tests | Collect success/failure matrix |

**Exit criteria:** Identify which geometries vertical MVP handles and which it fails on.

### Phase 2: Create MVP_NAT with Control Logic

**Goal:** Build `mvp_nat.py` plugin with:
- Vertical-first resolution selector
- MMO-aware fallback to heading if vertical fails
- Ceiling-aware fallback to heading if vertical fails
- VS-stop for PZ intrusion prevention
- Enhanced bouncing prevention for parallel tracks

| Step | Action |
|---|---|
| 2.1 | Implement `mvp_nat.py` skeleton extending MVP |
| 2.2 | Add vertical-first selector with ceiling/MMO fallback |
| 2.3 | Test on G9 (ceiling) and G8 (MMO) — verify fallback triggers |
| 2.4 | Add VS-stop logic |
| 2.5 | Test on G2 and G4 — verify VS-stop prevents vertical intrusion |
| 2.6 | Add anti-bouncing for parallel tracks |
| 2.7 | Test on G3 and G10 — verify no lateral drift |
| 2.8 | Run all G1–G10 with MVP_NAT |

**Exit criteria:** All unit test geometries resolve without LoS. Lateral deviation < 5 NM.

### Phase 3: Enhanced Recovery & Priority

**Goal:** Implement 2-criteria recovery, tune priority rules for NAT mix (cruising vs. climbing).

| Step | Action |
|---|---|
| 3.1 | Implement 2-criteria `resumenav()` override |
| 3.2 | Test on long-duration conflicts (G7 shallow crossing) |
| 3.3 | Experiment with FF2 vs FF3 vs LAY1 priority codes |
| 3.4 | Test with mixed cruising/climbing traffic (B2 stress scenarios) |

### Phase 4: Scale to NAT Separation (14 NM)

**Goal:** Change PZR from 5 NM to 14 NM and re-validate.

| Step | Action |
|---|---|
| 4.1 | Set `ZONER 14` in test config |
| 4.2 | Re-run all G1–G10 unit tests |
| 4.3 | Adjust `DTLOOK` (may need 600–900s for 14 NM at M0.85) |
| 4.4 | Tune `RFACH` / `RFACV` for 14 NM resolution margin |
| 4.5 | Run stress test B2 at 14 NM |

**Key calculation:** At M0.85 (~250 m/s or ~480 kt) head-on (closing ~960 kt), 14 NM takes ≈ 53 seconds to close. A 300s lookahead detects at ~80 NM. For same-direction (closing ~20 kt), 14 NM takes ≈ 25 min → need `DTLOOK ≥ 1500s` for early detection.

### Phase 5: Full-Day NAT Validation

**Goal:** Run 1000-flight scenario with MVP_NAT at 14 NM.

| Step | Action |
|---|---|
| 5.1 | Create NAT-border scenario from `optimised.scn` |
| 5.2 | Run with ASAS off → baseline conflict count |
| 5.3 | Run with MVP_NAT → measure LoS rate, fuel penalty, deviations |
| 5.4 | Analyse results with postprocessing pipeline |
| 5.5 | Iterate on parameters (DTLOOK, RFACH, priority code) |

### Phase 6: Intent-Based Enhancement (Future)

| Step | Action |
|---|---|
| 6.1 | Prototype intent-based conflict predictor using `traf.ap.route` |
| 6.2 | Compare intent-based vs. state-based detection accuracy |
| 6.3 | Integrate as optional enhancement in MVP_NAT |

---

## Summary: Immediate Next Actions

1. **Fix** `detection.py` line 10: `asas_pzr=1.0` → `asas_pzr=5.0`
2. **Create** `scenario/ASAS_NAT/` directory with unit test scenarios (G1–G10)
3. **Create** `test_ASAS_config.scn` with common ASAS setup
4. **Run** Phase 0 baseline detection tests
5. **Create** `bluesky/plugins/asas/mvp_nat.py` skeleton
6. **Enhance** `loggercompact.py` with ASAS diagnostic columns

---

## Appendix: ASAS Commands Quick Reference

| Command | Description | NAT Setting |
|---|---|---|
| `CDMETHOD STATEBASED` | Select state-based conflict detection | Always on |
| `RESO MVP` (or `MVP_NAT`) | Select resolution method | MVP_NAT when ready |
| `ZONER 5` / `ZONER 14` | Horizontal PZ radius [NM] | 5 → 14 |
| `ZONEDH 1000` | Vertical PZ half-height [ft] | 1000 |
| `DTLOOK 300` / `DTLOOK 1500` | Lookahead time [s] | 300 → 1500 |
| `RMETHH BOTH` | Enable heading + speed resolution | Depends on phase |
| `RMETHV ON` | Enable vertical resolution | Primary for NAT |
| `PRIORULES ON FF3` | Climbing/descending has priority | Recommended for NAT |
| `RFACH 1.05` | Horizontal resolution margin factor | 1.05 |
| `RFACV 1.05` | Vertical resolution margin factor | 1.05 |
| `NORESO acid` | Mark aircraft nobody avoids | For testing |
| `RESOOFF acid` | Mark aircraft that won't avoid | For testing |
| `STARTLOG` | Start compact logger | Always first |