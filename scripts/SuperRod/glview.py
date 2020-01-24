import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np

color_lib = {'C':(1,0,0,1),'O':(0,1,0,1),'Cu':(1,0,1,1)}
class GLViewWidget_cum(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

    def show_structure(self, xyz):
        self.setCameraPosition(distance=55, azimuth=-90)
        for each in xyz:
            e, x, y, z = each
            g = gl.GLGridItem()
            g.scale(2,2,1)
            self.addItem(g)

            md = gl.MeshData.sphere(rows=10, cols=20)

            m1 = gl.GLMeshItem(meshdata=md, smooth=True, color=color_lib[e], shader='shaded', glOptions='opaque')
            m1.translate(x, y, z)
            m1.scale(0.5, 0.5, 0.5)
            self.addItem(m1)

    def setup_view(self):
        self.setCameraPosition(distance=15, azimuth=-90)
        g = gl.GLGridItem()
        g.scale(2,2,1)
        self.addItem(g)

        md = gl.MeshData.sphere(rows=10, cols=20)
        x = np.linspace(-8, 8, 6)

        m1 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0), shader='balloon', glOptions='additive')
        m1.translate(x[0], 0, 0)
        m1.scale(1, 1, 2)
        self.addItem(m1)

        m2 = gl.GLMeshItem(meshdata=md, smooth=True, shader='normalColor', glOptions='opaque')
        m2.translate(x[1], 0, 0)
        m2.scale(1, 1, 2)
        self.addItem(m2)

        m3 = gl.GLMeshItem(meshdata=md, smooth=True, shader='viewNormalColor', glOptions='opaque')
        m3.translate(x[2], 0, 0)
        m3.scale(1, 1, 2)
        self.addItem(m3)

        m4 = gl.GLMeshItem(meshdata=md, smooth=True, shader='shaded', glOptions='opaque')
        m4.translate(x[3], 0, 0)
        m4.scale(1, 1, 2)
        self.addItem(m4)

        m5 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0, 1), shader='edgeHilight', glOptions='opaque')
        m5.translate(x[4], 0, 0)
        m5.scale(1, 1, 2)
        self.addItem(m5)

        m6 = gl.GLMeshItem(meshdata=md, smooth=True, color=(1, 0, 0, 1), shader='heightColor', glOptions='opaque')
        m6.translate(x[5], 0, 0)
        m6.scale(1, 1, 2)
        self.addItem(m6)