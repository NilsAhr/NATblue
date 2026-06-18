# MVP Resolution Method — Test Scenarios

These scenarios systematically test the Modified Voltage Potential (MVP) conflict resolution method under different configurations.

## Quick Reference: MVP Commands

| Command | Description |
|---------|-------------|
| `RESO MVP` | Activate MVP as resolution method |
| `ASAS ON` | Activate conflict detection |
| `RMETHH ON/BOTH` | Horizontal resolution: heading + speed (vertical OFF) |
| `RMETHH HDG` | Horizontal resolution: heading only |
| `RMETHH SPD` | Horizontal resolution: speed only |
| `RMETHH OFF` | Disable horizontal limitation |
| `RMETHV ON` or `RMETHV V/S` | Vertical resolution only (horizontal OFF) |
| `RMETHV OFF` | Disable vertical limitation |
| `RFACH <factor>` | Horizontal resolution margin (default 1.01) |
| `RFACV <factor>` | Vertical resolution margin (default 1.01) |
| `PRIORULES ON <code>` | Enable priority rules (FF1/FF2/FF3/LAY1/LAY2) |

## Test Scenarios

| # | File | Configuration | What it tests |
|---|------|--------------|---------------|
| 01 | `MVP_01_default_horiz.scn` | Default (SPD+HDG) | Baseline head-on horizontal resolution |
| 02 | `MVP_02_hdg_only.scn` | `RMETHH HDG` | Heading-only resolution |
| 03 | `MVP_03_spd_only.scn` | `RMETHH SPD` | Speed-only resolution |
| 04 | `MVP_04_vertical_only.scn` | `RMETHV ON` | Vertical-only resolution |
| 05 | `MVP_05_combined_horiz_vert.scn` | Both OFF | Combined horizontal + vertical |
| 06 | `MVP_06_crossing_90deg.scn` | Default | 90-degree crossing conflict |
| 07 | `MVP_07_overtaking.scn` | Default | Overtaking (same direction) conflict |
| 08 | `MVP_08_prio_FF1.scn` | `PRIORULES ON FF1` | No priority, cooperative |
| 09 | `MVP_09_prio_FF2.scn` | `PRIORULES ON FF2` | Cruising aircraft has priority |
| 10 | `MVP_10_multi_aircraft.scn` | Default | 4-way conflict (MVP limitation test) |
| 11 | `MVP_11_resolution_margin.scn` | `RFACH 1.5` | Larger resolution margins |

## How to run

In BlueSky, use:
```
PCALL scenario/MVP_tests/MVP_01_default_horiz.scn
```
Or open BlueSky and use `IC scenario/MVP_tests/MVP_01_default_horiz.scn`.

## Mode interaction matrix

| RMETHH | RMETHV | Channels used |
|--------|--------|---------------|
| (default/ON) | (default/OFF) | HDG + SPD |
| HDG | OFF | HDG only |
| SPD | OFF | SPD only |
| OFF | ON | V/S only |
| OFF | OFF | HDG + SPD + V/S |

> **Important**: RMETHH and RMETHV are mutually exclusive.  
> Setting RMETHH ON automatically sets RMETHV OFF and vice versa.  
> Setting both to OFF enables full 3D resolution.
