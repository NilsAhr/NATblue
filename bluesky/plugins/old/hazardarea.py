"""
BlueSky dynamic hazard area plugin.
This plugin can use an area definition to evacuate aircraft from that area.

Created on 2025 March
@author: ahre_ni, Nils Ahrenhold (TUD/DLR)

"""
import numpy as np
from bluesky import stack, traf, sim
from bluesky.tools import areafilter, geo
from bluesky.core import Entity, timed_function
import bluesky as bs

def init_plugin():
    dha = DynamicHazardArea()
    config = {
        'plugin_name': 'HAZARDAREA',
        'plugin_type': 'sim'
    }
    
    """
    stackfunctions = {
        'ACTIVATEDHA': [
            'ACTIVATEDHA name,method',
            'txt,txt',
            dha.activatedynamichazardarea,
            'Activate a dynamic hazard area with a specified evacuation method.'
        ],
        'DEACTIVATEDHA': [
            'DEACTIVATEDHA',
            '',
            dha.deactivate,
            'Deactivate the current dynamic hazard area.'
        ]
    }
    """
    
    return config #, stackfunctions

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
    """
    # Activates a DHA and sets the evacuation method
    def activatedynamichazardarea(self, name, method="SP"):
        if not areafilter.hasArea(name):
            return False, f"DHA {name} does not exist."
        
        self.active = True
        self.dha_name = name
        self.evacuation_method = method
        self.inside_dha.clear()
        self.outside_dha.clear()
        return True, f"Dynamic Hazard Area: {name} activated with method: {method}."
    
    # Deactivates the DHA and resets all aircraft states
    def deactivate(self):
        self.active = False
        self.dha_name = ""
        self.inside_dha.clear()
        self.outside_dha.clear()
        return True, "DHA deactivated."
    """
    def getPolygonEdges(self):
        dha_polygon = areafilter.basic_shapes[self.dha_name].coordinates
        # returning the polygon line segments based on the edges
        dha_edges = [(dha_polygon[i], dha_polygon[i + 1], dha_polygon[i + 2], dha_polygon[i + 3]) 
                     for i in range(0, len(dha_polygon) - 2, 2)]
        # Add the last segment to close the polygon
        dha_edges.append((dha_polygon[-2], dha_polygon[-1], dha_polygon[0], dha_polygon[1]))
        return dha_edges

    @timed_function(name='DHA', dt=1.0)
    #Checks aircraft positions and applies evacuation or holding.
    def update(self, dt):
        if not self.active:
            return
        
        lat, lon, alt = traf.lat, traf.lon, traf.alt
        in_dha = areafilter.checkInside(self.dha_name, lat, lon, alt)
        
        self.inside_dha = [id for id, inside in zip(traf.id, in_dha) if inside]
        self.outside_dha = [id for id, inside in zip(traf.id, in_dha) if not inside]
        
        dha_edges = self.getPolygonEdges() #should be handed to exevute evacuation

        for acid in self.inside_dha:
            self.execute_evacuation(acid, dha_edges)
        
        for acid in self.outside_dha:
            self.execute_holding(acid)

        self.active = False #only execute once the evacuation

    # Currently implements "Speed Up" (SU). Other methods will be added.  
    def execute_evacuation(self, acid, dha_edges):
        idx = traf.id2idx(acid)
        if self.evacuation_method == "SU":
            bs.scr.echo("Evacuation: speed up (SU)")
            stack.stack(f"{acid} spd {self.maxspeed}") #limits the speed to performance, see traffic "perf" based on perfbase.py
            #bs.traf.selspd[idx] = self.maxspeed #limits the speed to performance, see traffic "perf" based on perfbase.py
        elif self.evacuation_method == "SP":
            bs.scr.echo("Evacuation: shortest path (SP)")
            exit_heading = self.calculate_shortest_exit(acid, dha_edges)
            stack.stack(f"{acid} hdg {exit_heading}")
        elif self.evacuation_method == "SPSU":
            bs.scr.echo("Evacuation: shortest path and speed up (SPSU)")
            stack.stack(f"{acid} spd {self.maxspeed}")
            exit_heading = self.calculate_shortest_exit(acid, dha_edges)
            stack.stack(f"{acid} hdg {exit_heading}")
    
    
    # Assigns a holding pattern to aircraft outside the DHA.
    def execute_holding(self, acid):
        stack.stack(f"{acid} spd {self.minspeed}")
        #idx = traf.id2idx(acid)
        #bs.traf.selspd[idx] = self.minspeed

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
    
    @stack.command
    def activatedha(self, name: str, method: str):
        """ Activate a dynamic hazard area with a specified evacuation method. """
        if not areafilter.hasArea(name):
            return False, f"DHA {name} does not exist."
        
        self.active = True
        self.dha_name = name
        
        if not method:
            self.evacuation_method = "SP"
        self.evacuation_method = method
        self.inside_dha.clear()
        self.outside_dha.clear()
        return True, f"Dynamic Hazard Area: {name} activated with method: {method}."


    @stack.command
    def deactivatedha(self, name=str):
        """ Deactivate the current dynamic hazard area. """
        self.active = False
        self.dha_name = name
        self.inside_dha.clear()
        self.outside_dha.clear()
        return True, f"DHA: {name} deactivated."        

    
