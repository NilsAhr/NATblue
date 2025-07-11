""" Bird traffic gui plugin """
import numpy as np

import bluesky as bs
from bluesky import ui 
from bluesky.ui import palette
from bluesky import settings
import bluesky.ui.qtgl.glhelpers as glh
from bluesky.tools.aero import ft
from bluesky.ui.qtgl.guiclient import UPDATE_ALL

# Register settings defaults
settings.set_variable_defaults(text_size=13, bird_size=10)

palette.set_default_colours(
    bird=(255, 255, 0))

# Static
MAX_NBIRDS = 10000

### Initialization function of your plugin.
def init_plugin():
    config = {
        'plugin_name':     'BIRDGUI',
        'plugin_type':     'gui',
        }

    return config

# Bird traffic class
class BirdTraffic(ui.RenderObject, layer=100):
    # TODO: FIX labels
    def __init__(self, parent):
        super().__init__(parent=parent)
        
        self.bird_trk = glh.GLBuffer()
        self.bird_lat = glh.GLBuffer()
        self.bird_lon = glh.GLBuffer()
        self.bird_alt = glh.GLBuffer()
        self.bird_color = glh.GLBuffer()
        self.bird_lbl = glh.GLBuffer()
        self.bird_symbol = glh.VertexArrayObject(glh.gl.GL_TRIANGLE_FAN)
        self.birdlabels = glh.Text(settings.text_size, (8, 3))
        self.nbirds = 0
        self.show_lbl = True

        # subscribe to BIRDDATA stream only from active node
        bs.net.subscribe(b'BIRDDATA', actonly=True)
        
        # connect the stream to bird catcher
        bs.net.stream_received.connect(self.bird_catcher)

        # get stream of actnodedata changed to reset birds
        bs.net.actnodedata_changed.connect(self.bird_reset)

    def create(self):
        
        bird_size = settings.bird_size
        self.bird_trk.create(MAX_NBIRDS * 4, glh.GLBuffer.StreamDraw)
        self.bird_lat.create(MAX_NBIRDS * 4, glh.GLBuffer.StreamDraw)
        self.bird_lon.create(MAX_NBIRDS * 4, glh.GLBuffer.StreamDraw)
        self.bird_alt.create(MAX_NBIRDS * 4, glh.GLBuffer.StreamDraw)
        self.bird_color.create(MAX_NBIRDS * 4, glh.GLBuffer.StreamDraw)
        self.bird_lbl.create(MAX_NBIRDS * 24, glh.GLBuffer.StreamDraw)


        # ------- Bird triangle fan -------------------------
        birdvertices = np.array([(0.0, 0.0),                     # 1
                            (0.7 * bird_size, -0.1 * bird_size),   # wing, 2 
                            (0.7 * bird_size, 0.1 * bird_size),    # wing, 3
                            (0.25 * bird_size, 0.1 * bird_size),    # trans, 4
                            (0.1 * bird_size, 0.1 * bird_size),    # trans 5
                            (0.05 * bird_size, 0.5 * bird_size),    #head. 6
                            (-0.05 * bird_size, 0.5 * bird_size),    #head. 7
                            (-0.1 * bird_size, 0.1 * bird_size),    # trans 8
                            (-0.25 * bird_size, 0.1 * bird_size),    # trans, 9
                            (-0.7 * bird_size, 0.1 * bird_size),    # wing, 10
                            (-0.7 * bird_size, -0.1 * bird_size),   # wing, 11
                            (-0.25 * bird_size, -0.1 * bird_size),   # trans, 12 
                            (0, -0.5 * bird_size),   # tail, 13
                            (0.25 * bird_size, -0.1 * bird_size),   # trans, 14 
                            (0.7 * bird_size, -0.1 * bird_size)],
                        dtype=np.float32)

        self.bird_symbol.create(vertex=birdvertices)

        self.bird_symbol.set_attribs(lat=self.bird_lat, lon=self.bird_lon, color=self.bird_color,
                                   orientation=self.bird_trk, instance_divisor=1)

        self.birdlabels.create(self.bird_lbl, self.bird_lat, self.bird_lon, self.bird_color,
                             (bird_size, -0.5 * bird_size), instanced=True)

    def draw(self):
        #print ('GL - bird drawing')
        if self.nbirds:
            self.bird_symbol.draw(n_instances=self.nbirds)

        # TODO: fix label
        # if self.show_lbl:
        #     self.birdlabels.draw(n_instances=self.nbirds)

    def update_bird_data(self, data):
        #print ('GL - bird updating', data['id'])
        
        # get bird data
        bird_id = data['id']
       # bird_type = data['type']
        bird_lat = data['lat']
        bird_lon = data['lon']
        bird_trk = data['trk']
        bird_alt = data['alt']
       # bird_vs = data['vs']
       # bird_hs = data['hs']

        # update buffers
        self.nbirds = len(bird_lat)
        self.bird_lat.update(np.array(bird_lat, dtype=np.float32))
        self.bird_lon.update(np.array(bird_lon, dtype=np.float32))
        self.bird_trk.update(np.array(bird_trk, dtype=np.float32))
        self.bird_alt.update(np.array(bird_alt, dtype=np.float32))

        # colors
        rawlabel = ''
        color = np.empty(
            (min(self.nbirds, MAX_NBIRDS), 4), dtype=np.uint8)
        rgb_bird = palette.bird

        zdata = zip(data['id'], data['alt'])
        for i, (id, alt) in enumerate(zdata):
            if i >= MAX_NBIRDS:
                break

            # Make label
            if self.show_lbl:
                #print ("in gl, birdie id is", id, type(id), len(str(id)))
                #print ("spooky abbreviation", '%-8s')
                rawlabel += '%-8s' % str(id)[:8]
                rawlabel += '%-5d' % int(alt / ft + 0.5)

            color[i, :] = tuple(rgb_bird) + (255,)
        
        # update bird label
        self.bird_color.update(color)
        self.bird_lbl.update(np.array(rawlabel.encode('utf8'), dtype=np.string_))

    def bird_catcher(self, name, data, sender_id):
        """receive stream from bluesky sim.
        """
        #print ('GL - catching birdies')
        # update bird data if stream is bird data
        if name == b'BIRDDATA':
            self.update_bird_data(data)

    def bird_reset(self, nodeid, nodedata, changed_elems):
        """Receive signal from actnodedata_changed.
        When changed_elems is equal to UPDATE_ALL from gui client 
        it signals a reset or a new node.
        """
        #print('GL - bird reset is called')
        # reset birds if changing to a new active node
        if changed_elems == UPDATE_ALL:
            
            data = dict()

            data['id']         =[]
            data['type']       = []
            data['lat']        = []
            data['lon']        = []
            data['alt']        = []
            data['trk']        = []
            data['vs']         = []
            data['hs']         = []

            self.update_bird_data(data)
        