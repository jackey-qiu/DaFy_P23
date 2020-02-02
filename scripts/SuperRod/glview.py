import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from pyqtgraph.Qt import QtGui

color_lib = {'C':(1,0,0,1),'O':(0,1,0,1),'Cu':(1,0,1,1)}
class GLViewWidget_cum(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

    def show_structure(self, xyz):
        # self.setCameraPosition(distance=55, azimuth=-90)
        self.setCameraPosition(distance=55, azimuth=0)
        i=0
        if len(self.items)==0:
            xgrid = gl.GLGridItem()
            ygrid = gl.GLGridItem()
            zgrid = gl.GLGridItem()

            xgrid.setSize(5,5,5)
            ygrid.setSize(5,5,5)
            zgrid.setSize(5,5,5)
 
            self.addItem(xgrid)
            self.addItem(ygrid)
            self.addItem(zgrid)
 
            # Rotate x and y grids to face the correct direction
            xgrid.rotate(90, 0, 1, 0)
            ygrid.rotate(90, 1, 0, 0)
 
            # Scale grids to the appropriate dimensions
            xgrid.scale(3.15, 3.15, 3.15)
            ygrid.scale(3.15, 3.15, 3.15)
            zgrid.scale(3.15, 3.15, 3.15)
            """
            g = gl.GLGridItem(size = QtGui.QVector3D(30,30,1))
            g.scale(1,1,1)
            # g.setSpacing(3.54,3.54,3.54)
            g.rotate(90,1,0,0)
            self.addItem(g)
            g1 = gl.GLGridItem(size = QtGui.QVector3D(30,30,1))
            g1.scale(1,1,1)
            g1.setSpacing(3.54,3.54,3.54)
            #g1.rotate(90,1,0,0)
            self.addItem(g1)
            """
            for each in xyz:
                e, x, y, z = each
                md = gl.MeshData.sphere(rows=10, cols=20)
                m1 = gl.GLMeshItem(meshdata=md, smooth=True, color=color_lib[e], shader='shaded', glOptions='opaque')
                # print(dir(m1.metaObject()))
                m1.translate(x, y, z)
                m1.scale(0.5, 0.5, 0.5)
                self.addItem(m1)
        else:
            for each in xyz:
                _,x,y,z = xyz[i]
                #first item is grid net
                self.items[i+1+1+1].resetTransform()
                self.items[i+1+1+1].translate(x,y,z)
                i += 1

    def update_structure(self, xyz):
        for i in range(len(xyz)):
            _,x,y,z = xyz[i]
            #first item is grid net
            self.items[i+1+1+1].resetTransform()
            self.items[i+1+1+1].translate(x,y,z)
            self.items[i+1+1+1].scale(0.5, 0.5, 0.5)

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