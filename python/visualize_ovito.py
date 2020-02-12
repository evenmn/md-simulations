import math
from ovito.io import import_file
from ovito.vis import Viewport

def visualize(infile, outfile='out.png'):
    '''
    This function takes a xyz-file as an argument and visualizes
    the particles using Ovito.
    '''
    pipeline = import_file(infile)      # Create a pipeline object
    pipeline.add_to_scene()             # Add pipeline to three-dimensional scene
    vp = Viewport()                     # Object defining viewpoint from which the scene is seen
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (-250, -250, 300)
    vp.camera_dir = (3, 3, -3)
    vp.fov = math.radians(60.0)
    vp.render_image(size=(800,600), filename=outfile, background=(0,0,0), frame=8)
    pipeline.remove_from_scene()
