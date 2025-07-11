"""
BlueSky dynamic hazard area plugin.
This plugin can use an area definition to evacuate aircraft from that area.

Created on 2025 March
@author: ahre_ni, Nils Ahrenhold (TUD/DLR)

"""
import numpy as np
from bluesky import core, stack, traf, sim
from bluesky.tools import areafilter, geo
from bluesky.core import Entity, timed_function
import bluesky as bs


def init_plugin():
    global dynamichazardarea
    dynamichazardarea = DynamicHazardArea()

    config = {
        'plugin_name': 'DYNAMICHAZARDAREA',
        'plugin_type': 'sim'
    }
        
    return config

class DynamicHazardArea(Entity):

    def __init__(self):
        super().__init__()

        self.active = False
        self.dha_name = ""
        self.inside_dha = []  # Aircraft inside DHA
        self.outside_dha = []  # Aircraft outside DHA
        self.evacuation_method = "SU"  # Default: Speed Up
        self.maxspeed = 800
        self.minspeed = 0
        self.min_dist = float("inf")
    
    def getPolygonEdges(self):
        dha_polygon = areafilter.basic_shapes[self.dha_name].coordinates
        # returning the polygon line segments based on the edges
        dha_edges = [(dha_polygon[i], dha_polygon[i + 1], dha_polygon[i + 2], dha_polygon[i + 3]) 
                     for i in range(0, len(dha_polygon) - 2, 2)]
        # Add the last segment to close the polygon
        dha_edges.append((dha_polygon[-2], dha_polygon[-1], dha_polygon[0], dha_polygon[1]))
        return dha_edges

    # Currently implements "Speed Up" (SU). Other methods will be added.  
    def execute_evacuation(self, acid, dha_edges, method):
        if method == "SU":
            bs.scr.echo("Evacuation: speed up (SU)")
            stack.stack(f"{acid} spd {self.maxspeed}") #limits the speed to performance, see traffic "perf" based on perfbase.py
        elif method == "SP":
            bs.scr.echo("Evacuation: shortest path (SP)")
            exit_heading = self.calculate_shortest_exit(acid, dha_edges)
            stack.stack(f"{acid} hdg {exit_heading}")
        elif method == "SPSU":
            bs.scr.echo("Evacuation: shortest path and speed up (SPSU)")
            stack.stack(f"{acid} spd {self.maxspeed}")
            exit_heading = self.calculate_shortest_exit(acid, dha_edges)
            stack.stack(f"{acid} hdg {exit_heading}")
        else :
            bs.scr.echo("Evacuation: without change (NO)")
            
    
    # Assigns a holding pattern to aircraft outside the DHA.
    def execute_holding(self, acid):
        stack.stack(f"{acid} spd {self.minspeed}")

    def calculate_shortest_exit(self, acid, dha_edges):
        idx = traf.id2idx(acid)
        ac_lat, ac_lon = traf.lat[idx], traf.lon[idx]   
        min_dist = self.min_dist
        best_heading = traf.hdg[idx]

        for lat1, lon1, lat2, lon2 in dha_edges:
            proj_lat, proj_lon = geo.project_point_on_line(ac_lat, ac_lon, lat1, lon1, lat2, lon2)
            qdr, dist = geo.qdrdist(ac_lat, ac_lon, proj_lat, proj_lon)
            if dist < min_dist:
                min_dist = dist
                best_heading = qdr
                best_heading = (best_heading + 360) % 360
        bs.scr.echo(f"best heading for {acid} is : {best_heading}")
        
        return best_heading
    
    @core.timed_function(name='dynamichazardarea', dt=1)
    def update(self):
        #Checks if all aircraft left DHA
        n = traf.ntraf
        if n == 0:
            stack.stack("hold")
            stack.stack("close")

    @stack.command
    def activatedha(self, name: str, method: str = None):
        """ Activate a dynamic hazard area with a specified evacuation method. """
        if not areafilter.hasArea(name):
            return False, f"DHA {name} does not exist."
        
        self.active = True
        self.dha_name = name       
        self.evacuation_method = method if method else "SP"  # Default to SP if empty
        self.inside_dha.clear()
        self.outside_dha.clear()
        # retrieve aircraft information
        lat, lon, alt = traf.lat, traf.lon, traf.alt
        in_dha = areafilter.checkInside(self.dha_name, lat, lon, alt)
        
        self.inside_dha = [id for id, inside in zip(traf.id, in_dha) if inside]
        self.outside_dha = [id for id, inside in zip(traf.id, in_dha) if not inside]
        bs.scr.echo(f"inside is: {self.inside_dha}")
        dha_edges = self.getPolygonEdges() #should be handed to exevute evacuation

        for acid in self.inside_dha:
            bs.scr.echo("entering loop evac")
            self.execute_evacuation(acid, dha_edges, method)
            bs.scr.echo(f"execute evac for {acid} with: {method}")
        
        for acid in self.outside_dha:
            #self.execute_holding(acid)
            stack.stack(f"{acid} del")

        return True, f"Dynamic Hazard Area: {name} activated with method: {method}."

        
    @stack.command
    def deactivatedha(self, name=str):
        """ Deactivate the current dynamic hazard area. """
        self.active = False
        self.dha_name = name
        self.inside_dha.clear()
        self.outside_dha.clear()
        return True, f"DHA: {name} deactivated."

    
        
        