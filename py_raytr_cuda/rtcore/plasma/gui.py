#!/usr/bin/python

import raytrace as rt

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


from pylab import *
import numpy as nu
from matplotlib.colors import LogNorm

import os

#from enthought.mayavi.mlab import *

import pygtk
pygtk.require('2.0')
import gtk

class MainWindow:

    # This is a callback function. The data arguments are ignored
    # in this example. More on callbacks below.
    def trace(self, widget, data=None):
        grid = (self.grid_spinner_width.get_value(),
                self.grid_spinner_height.get_value())
        rect = (self.implane_spinner_xMin.get_value(),
                self.implane_spinner_yMin.get_value(),
                self.implane_spinner_xMax.get_value(),
                self.implane_spinner_yMax.get_value())
        obs = self.observer_field.get_text().split(',')
        sphere = float(self.sphere_field.get_text())
        freq = float(self.freq_field.get_text())
        if(self.mode_radio_basic.get_active()):
            mode = 'basic'
            self.graph_button_tbr.set_sensitive(False)
            self.graph_button_tbriv.set_sensitive(False)
        else:
            if(self.mode_radio_tbr.get_active()):
                mode = 'Tbr'
                self.graph_button_tbr.set_sensitive(True)
                self.graph_button_tbriv.set_sensitive(False)
            else:
                mode = 'TbrIV'
                self.graph_button_tbriv.set_sensitive(True)
                self.graph_button_tbr.set_sensitive(False)
                
                
        cnu = float(self.cnu_field.get_text())
        msun = float(self.msun_field.get_text())
        
        self.trace = rt.implane(grid,
                       rect,
                       (float(obs[0]),float(obs[1]),float(obs[2])),
                       sphere,
                       freq,
                       mode,
                       cnu,
                       msun)
        self.trace.set_plfunc('plasma_parameters.c')
        if(self.streamerSet!=None):
            for streamer in self.streamerSet:
                self.trace.make_streamer(streamer[0]*180/math.pi,
                                         streamer[1]*180/math.pi,
                                         orientation=-streamer[2]+math.pi/2)
            self.streamerSet.clear()

        if(len(self.selected_set)!=0):
            trajPoints = empty((len(self.selected_set),2))
            i = 0
            for selected in self.selected_set:
                trajPoints[i] = selected
                i = i + 1
            self.traj = self.trace.trace(1500,trajPoints)
            self.graph_button_traj.set_sensitive(True)
            
            

            

        else:
            self.trace.trace(1500)
            self.graph_button_traj.set_sensitive(False)
        

    def tbr(self, widget, data=None):
        fig = figure()
        imshow(self.trace.tbr)
        colorbar()
        #fig.suptitle('Tbr, frequency = '+self.freq_field.getText())
        show()

    def tbriv(self, widget, data=None):
        fig = figure()
        imshow(self.trace.tbriv[:,:,0])
        colorbar()
        #fig.suptitle('TbrI, frequency = '+self.freq_field.getText())
        show()

        fig = figure()
        imshow(self.trace.tbriv[:,:,1])
        colorbar()
        #fig.suptitle('TbrV, frequency = '+self.freq_field.getText())
        show()

    def trajectory(self, widget, data=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')


        
        for trajectory in self.traj[:]:
            ax.plot(trajectory[:,0], trajectory[:,1], trajectory[:,2])

        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 20)

        x = 1 * np.outer(np.cos(u), np.sin(v))
        y = 1 * np.outer(np.sin(u), np.sin(v))
        z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='y')

        #fig.suptitle('Trajectories, frequency = '+self.freq_field.getText())

        plt.show()

    def clear(self, widget, data=None):
        self.selected_set.clear()
        self.sun_drawing.queue_draw()

    def selectAll(self, widget, data=None):
        for i in range(int(self.grid_spinner_width.get_value())):
            for j in range(int(self.grid_spinner_height.get_value())):
                self.selected_set.add((i,j))
        self.sun_drawing.queue_draw()

    def streamers(self, widget, data=None):
        self.streamerSet = StreamerWindow(self.window, self.streamerSet).getStreamers()
        

    def sunExpose(self, widget, data=None):
        style = self.sun_drawing.get_style()
        gc = style.fg_gc[gtk.STATE_NORMAL]

        gc.set_rgb_fg_color(gtk.gdk.Color('#ff0'))

        width = self.sun_drawing.size_request()[0]
        height = self.sun_drawing.size_request()[1]

        xLength = float(self.implane_spinner_xMax.get_value()) - float(self.implane_spinner_xMin.get_value())
        yLength = float(self.implane_spinner_yMax.get_value()) - float(self.implane_spinner_yMin.get_value())

        xPos = -float(self.implane_spinner_xMin.get_value())/xLength*width
        yPos = -float(self.implane_spinner_yMin.get_value())/yLength*height

        xSize = 2.0/xLength*width
        ySize = 2.0/yLength*height

        self.sun_drawing.window.draw_arc(gc,True, int(xPos-xSize/2.0),int(yPos-ySize/2.0),int(xSize),int(ySize),0,360*64)

        xSpacing = (width-1)/float(self.grid_spinner_width.get_value())
        ySpacing = (height-1)/float(self.grid_spinner_height.get_value())

        gc.set_rgb_fg_color(gtk.gdk.Color('#44f'))

        x = int(self.highlightBox[0]*xSpacing)
        y = int(self.highlightBox[1]*ySpacing)
        w = int(xSpacing+1)
        h = int(ySpacing+1)

        self.sun_drawing.window.draw_rectangle(gc,True,x,y,w,h)         

        gc.set_rgb_fg_color(gtk.gdk.Color('#f44'))

        for box in self.selected_set:

            x = int(box[0]*xSpacing)
            y = int(box[1]*ySpacing)
            w = int(xSpacing+1)
            h = int(ySpacing+1)

            self.sun_drawing.window.draw_rectangle(gc,True,x,y,w,h)      

        gc.set_rgb_fg_color(gtk.gdk.Color('#000'))

        for i in arange(0,width,xSpacing):
            self.sun_drawing.window.draw_line(gc,int(i),0,int(i),int(height))
            
        for i in arange(0,height,ySpacing):
            self.sun_drawing.window.draw_line(gc,0,int(i),int(width),int(i))

    def mouseMotion(self, widget, data=None):

        width = self.sun_drawing.size_request()[0]
        height = self.sun_drawing.size_request()[1]

        xSpacing = (width-1)/float(self.grid_spinner_width.get_value())
        ySpacing = (height-1)/float(self.grid_spinner_height.get_value())
        
        oldHighlightBox = self.highlightBox
        self.highlightBox = (int(self.sun_drawing.get_pointer()[0]/xSpacing),
                             int(self.sun_drawing.get_pointer()[1]/ySpacing))
        if(oldHighlightBox == self.highlightBox):
            return True

        if(self.mouseHeld==1):
            self.selected_set.add(self.highlightBox)
        if(self.mouseHeld==2):
            if(self.highlightBox in self.selected_set):
                self.selected_set.remove(self.highlightBox)
        
        self.sun_drawing.queue_draw()
        
    def mouseExit(self, widget, data=None):
        self.highlightBox = (-1,-1)
        self.sun_drawing.queue_draw()

    def mouseClicked(self, widget, data=None):
        if(self.highlightBox != (-1,-1)):
            if(data.button==1):
                self.selected_set.add(self.highlightBox)
                self.sun_drawing.queue_draw()
                self.mouseHeld = 1
            if(data.button==2 or data.button==3):
                self.mouseHeld = 2
                if(self.highlightBox in self.selected_set):
                    self.selected_set.remove(self.highlightBox)
                    self.sun_drawing.queue_draw()

    def mouseReleased(self, widget, data=None):
        self.mouseHeld = 0

    def paramsChange(self, widget, data=None):
        self.selected_set.clear()
        self.sun_drawing.queue_draw()

    def save(self, widget, data=None):
        saver = gtk.FileChooserDialog(title='Save As', parent=self.window,
                              action=gtk.FILE_CHOOSER_ACTION_SAVE,
                              buttons=(gtk.STOCK_SAVE, gtk.RESPONSE_OK,
                                       gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))
        response = saver.run()
        if response == gtk.RESPONSE_OK:
            self.trace.save_fits(saver.get_filename())
        saver.destroy()

    def makeBatch(self, widget, data=None):
        dialog = gtk.Dialog(title='Make Batch Job',
                            parent=self.window,
                            buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK,
                                     gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))

        freq_box = gtk.HBox(False,10)

        freq_label_start = gtk.Label("Frequency Start:")
        freq_label_start.show()

        freq_box.pack_start(freq_label_start)
        
        freq_field_start = gtk.Entry(100)
        freq_field_start.set_text(self.freq_field.get_text())
        freq_field_start.show()
        
        freq_box.pack_start(freq_field_start)

        freq_label_end = gtk.Label("Frequency End:")
        freq_label_end.show()

        freq_box.pack_start(freq_label_end)
        
        freq_field_end = gtk.Entry(100)
        freq_field_end.set_text(self.freq_field.get_text())
        freq_field_end.show()
        
        
        freq_box.pack_start(freq_field_end)

        freq_box.show()

        freq_label_step = gtk.Label("Number of Frequency:")
        freq_label_step.show()

        freq_box.pack_start(freq_label_step)
        
        freq_field_step = gtk.Entry(100)
        freq_field_step.set_text('1')
        freq_field_step.show()
        
        
        freq_box.pack_start(freq_field_step)

        freq_box.show()
        
        dialog.vbox.pack_start(freq_box)



        theta_box = gtk.HBox(False,10)

        theta_label_start = gtk.Label("Theta Start:")
        theta_label_start.show()

        theta_box.pack_start(theta_label_start)
        
        theta_field_start = gtk.Entry(100)
        theta_field_start.set_text('0')
        theta_field_start.show()
        
        theta_box.pack_start(theta_field_start)

        theta_label_end = gtk.Label("Theta End:")
        theta_label_end.show()

        theta_box.pack_start(theta_label_end)
        
        theta_field_end = gtk.Entry(100)
        theta_field_end.set_text('90')
        theta_field_end.show()
        
        
        theta_box.pack_start(theta_field_end)

        theta_box.show()

        theta_label_step = gtk.Label("Number of Angles:")
        theta_label_step.show()

        theta_box.pack_start(theta_label_step)
        
        theta_field_step = gtk.Entry(100)
        theta_field_step.set_text('4')
        theta_field_step.show()
        
        
        theta_box.pack_start(theta_field_step)

        theta_box.show()
        
        dialog.vbox.pack_start(theta_box)

        
        phi_box = gtk.HBox(False,10)

        phi_label_start = gtk.Label("Phi Start:")
        phi_label_start.show()

        phi_box.pack_start(phi_label_start)
        
        phi_field_start = gtk.Entry(100)
        phi_field_start.set_text('0')
        phi_field_start.show()
        
        phi_box.pack_start(phi_field_start)

        phi_label_end = gtk.Label("Phi End:")
        phi_label_end.show()

        phi_box.pack_start(phi_label_end)
        
        phi_field_end = gtk.Entry(100)
        phi_field_end.set_text('90')
        phi_field_end.show()
        
        
        phi_box.pack_start(phi_field_end)

        phi_box.show()

        phi_label_step = gtk.Label("Number of Angles:")
        phi_label_step.show()

        phi_box.pack_start(phi_label_step)
        
        phi_field_step = gtk.Entry(100)
        phi_field_step.set_text('4')
        phi_field_step.show()
        
        
        phi_box.pack_start(phi_field_step)

        phi_box.show()
        
        dialog.vbox.pack_start(phi_box)

        response = dialog.run()

        grid = (self.grid_spinner_width.get_value(),
                self.grid_spinner_height.get_value())
        rect = (self.implane_spinner_xMin.get_value(),
                self.implane_spinner_yMin.get_value(),
                self.implane_spinner_xMax.get_value(),
                self.implane_spinner_yMax.get_value())
        obs = self.observer_field.get_text().split(',')
        sphere = float(self.sphere_field.get_text())
        if(self.mode_radio_basic.get_active()):
            mode = 'basic'
        else:
            if(self.mode_radio_tbr.get_active()):
                mode = 'Tbr'
            else:
                mode = 'TbrIV'
                
                
        cnu = float(self.cnu_field.get_text())
        msun = float(self.msun_field.get_text())
  
        
        if response == gtk.RESPONSE_OK:
            f = open('tmp.batch.py','w')
            f.write('import raytrace as rt\nimport numpy as np\nimport os\n')



            f.write('for theta in np.linspace('+\
                    theta_field_start.get_text()+\
                    ','+theta_field_end.get_text()+','+\
                    theta_field_step.get_text()+'):\n')

            f.write('\tfor phi in np.linspace('+\
                    phi_field_start.get_text()+\
                    ','+phi_field_end.get_text()+','+\
                    phi_field_step.get_text()+'):\n')

            f.write('\t\tos.mkdir(\'data/theta=\'+str(theta)+\'_phi=\'+str(phi))\n')
            
            f.write('\t\tfor i in np.linspace('+\
                         freq_field_start.get_text()+\
                         ','+freq_field_end.get_text()+','+\
                         freq_field_step.get_text()+'):\n')
            f.write('\t\t\ta = rt.implane('+str(grid)+','+str(rect)+','+\
                    str(obs)+','+str(sphere)+',i,\''+str(mode)+'\','+\
                    str(cnu)+','+str(msun)+')\n')

            f.write('\t\t\ta.make_streamer(theta,phi)\n')
            
            f.write('\t\t\ta.trace(1500)\n')
            f.write('\t\t\ta.save_fits(\'data/theta=\'+str(theta)+\'_phi=\'+str(phi)+\'/\'+str(i/1e6)+\'MHz\'+\'.fits\')\n')
            os.system('python tmp.batch.py')
            dialog.destroy()
        if response == gtk.RESPONSE_CANCEL:
            dialog.destroy()
        
                
    def delete_event(self, widget, event, data=None):
        # If you return FALSE in the "delete_event" signal handler,
        # GTK will emit the "destroy" signal. Returning TRUE means
        # you don't want the window to be destroyed.
        # This is useful for popping up 'are you sure you want to quit?'
        # type dialogs.
        print "delete event occurred"

        # Change FALSE to TRUE and the main window will not be destroyed
        # with a "delete_event".
        return False

    def destroy(self, widget, data=None):
        print "destroy signal occurred"
        gtk.main_quit()

    def __init__(self):
        # create a new window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        
        # When the window is given the "delete_event" signal (this is given
        # by the window manager, usually by the "close" option, or on the
        # titlebar), we ask it to call the delete_event () function
        # as defined above. The data passed to the callback
        # function is NULL and is ignored in the callback function.
        self.window.connect("delete_event", self.delete_event)
    
        # Here we connect the "destroy" event to a signal handler.  
        # This event occurs when we call gtk_widget_destroy() on the window,
        # or if we return FALSE in the "delete_event" callback.
        self.window.connect("destroy", self.destroy)
    
        # Sets the border width of the window.
        self.window.set_border_width(10)
        self.window.set_destroy_with_parent(False)

        self.highlightBox = (-1,-1)
        self.selected_set = set()
        self.mouseHeld = 0
        self.traj = None
        self.streamerSet = None

        self.big_box = gtk.VBox(False,10)



#makes the spinbuttons for the implane coordinates

        # 1st row of implane coordinates (label)

        self.implane_box_1 = gtk.HBox(False,10)
        self.implane_box_2 = gtk.HBox(False,10)
        self.implane_box_3 = gtk.HBox(False,10)
        
        self.implane_label = gtk.Label("Implane Coordinates:")
        self.implane_label.show()
        self.implane_box_1.pack_start(self.implane_label)
        self.implane_box_1.show()
        self.big_box.pack_start(self.implane_box_1)

        # 2nd row of implane coordinates (x)
        
        self.implane_label_xMin = gtk.Label("xMin coordinate")
        self.implane_label_xMin.show()

        self.implane_box_2.pack_start(self.implane_label_xMin)
        
        self.implane_spinner_xMin = gtk.SpinButton(gtk.Adjustment(-2,-100,100,.1,1),.1,2)
        self.implane_spinner_xMin.show()
        self.implane_spinner_xMin.connect("value-changed",self.paramsChange,None)
        
        self.implane_box_2.pack_start(self.implane_spinner_xMin)

        self.implane_label_xMax = gtk.Label("xMax coordinate")
        self.implane_label_xMax.show()

        self.implane_box_2.pack_start(self.implane_label_xMax)
        
        self.implane_spinner_xMax = gtk.SpinButton(gtk.Adjustment(2,-100,100,.1,1),.1,2)
        self.implane_spinner_xMax.show()
        self.implane_spinner_xMax.connect("value-changed",self.paramsChange)
        
        self.implane_box_2.pack_start(self.implane_spinner_xMax)

        self.implane_box_2.show()

        self.big_box.pack_start(self.implane_box_2)

        # 3rd row of implane coordinates (y)

        self.implane_label_yMin = gtk.Label("yMin coordinate")
        self.implane_label_yMin.show()

        self.implane_box_3.pack_start(self.implane_label_yMin)
        
        self.implane_spinner_yMin = gtk.SpinButton(gtk.Adjustment(-2,-100,100,.1,1),.1,2)
        self.implane_spinner_yMin.show()
        self.implane_spinner_yMin.connect("value-changed",self.paramsChange)
        
        self.implane_box_3.pack_start(self.implane_spinner_yMin)

        self.implane_label_yMax = gtk.Label("yMax coordinate")
        self.implane_label_yMax.show()

        self.implane_box_3.pack_start(self.implane_label_yMax)
        
        self.implane_spinner_yMax = gtk.SpinButton(gtk.Adjustment(2,-100,100,.1,1),.1,2)
        self.implane_spinner_yMax.show()
        self.implane_spinner_yMax.connect("value-changed",self.paramsChange)
        
        self.implane_box_3.pack_start(self.implane_spinner_yMax)

        self.implane_box_3.show()

        self.big_box.pack_start(self.implane_box_3)

#makes the grid width and height spiners

        # 1st row of grid size (label)

        self.grid_box_1 = gtk.HBox(False,10)
        self.grid_box_2 = gtk.HBox(False,10)
        
        self.grid_label = gtk.Label("Grid Size:")
        self.grid_label.show()
        self.grid_box_1.pack_start(self.grid_label)
        self.grid_box_1.show()
        self.big_box.pack_start(self.grid_box_1)

        # 2nd row of grid coordinates (x)
        
        self.grid_label_width = gtk.Label("Grid Width")
        self.grid_label_width.show()

        self.grid_box_2.pack_start(self.grid_label_width)
        
        self.grid_spinner_width = gtk.SpinButton(gtk.Adjustment(40,0,1000,1,1),1,0)
        self.grid_spinner_width.show()
        self.grid_spinner_width.connect("value-changed",self.paramsChange)
        
        self.grid_box_2.pack_start(self.grid_spinner_width)

        self.grid_label_height = gtk.Label("Grid Height")
        self.grid_label_height.show()

        self.grid_box_2.pack_start(self.grid_label_height)
        
        self.grid_spinner_height = gtk.SpinButton(gtk.Adjustment(40,0,1000,1,1),1,0)
        self.grid_spinner_height.show()
        self.grid_spinner_height.connect("value-changed",self.paramsChange)
        
        self.grid_box_2.pack_start(self.grid_spinner_height)

        self.grid_box_2.show()

        self.big_box.pack_start(self.grid_box_2)

#makes the observer location text box

        self.observer_box = gtk.HBox(False,10)

        self.observer_label = gtk.Label("Observer Location:")
        self.observer_label.show()

        self.observer_box.pack_start(self.observer_label)
        
        self.observer_field = gtk.Entry(100)
        self.observer_field.set_text("215.0,0,0")
        self.observer_field.show()
        
        self.observer_box.pack_start(self.observer_field)

        self.observer_box.show()
        
        self.big_box.pack_start(self.observer_box)

#makes the sphere text box

        self.sphere_box = gtk.HBox(False,10)

        self.sphere_label = gtk.Label("Sphere Location:")
        self.sphere_label.show()

        self.sphere_box.pack_start(self.sphere_label)
        
        self.sphere_field = gtk.Entry(100)
        self.sphere_field.set_text("25.0")
        self.sphere_field.show()
        
        self.sphere_box.pack_start(self.sphere_field)

        self.sphere_box.show()
        
        self.big_box.pack_start(self.sphere_box)

#makes the frequency text box

        self.freq_box = gtk.HBox(False,10)

        self.freq_label = gtk.Label("Frequency:")
        self.freq_label.show()

        self.freq_box.pack_start(self.freq_label)
        
        self.freq_field = gtk.Entry(100)
        self.freq_field.set_text("100e6")
        self.freq_field.show()
        
        self.freq_box.pack_start(self.freq_field)

        self.freq_box.show()
        
        self.big_box.pack_start(self.freq_box)

#makes the cnu text box

        self.cnu_box = gtk.HBox(False,10)

        self.cnu_label = gtk.Label("cnu:")
        self.cnu_label.show()

        self.cnu_box.pack_start(self.cnu_label)
        
        self.cnu_field = gtk.Entry(100)
        self.cnu_field.set_text("3.0")
        self.cnu_field.show()
        
        self.cnu_box.pack_start(self.cnu_field)

        self.cnu_box.show()
        
        self.big_box.pack_start(self.cnu_box)

#makes the msun text box

        self.msun_box = gtk.HBox(False,10)

        self.msun_label = gtk.Label("msun:")
        self.msun_label.show()

        self.msun_box.pack_start(self.msun_label)
        
        self.msun_field = gtk.Entry(100)
        self.msun_field.set_text("1.0")
        self.msun_field.show()
        
        self.msun_box.pack_start(self.msun_field)

        self.msun_box.show()
        
        self.big_box.pack_start(self.msun_box)
        
        #Makes the mode drop down

        self.mode_box = gtk.HBox(False,10)

        self.mode_label = gtk.Label("Mode:")
        self.mode_label.show()

        self.mode_box.pack_start(self.mode_label)
        
        self.mode_radio_basic = gtk.RadioButton(None,'basic')
        self.mode_radio_tbr = gtk.RadioButton(self.mode_radio_basic,'tbr')
        self.mode_radio_tbriv = gtk.RadioButton(self.mode_radio_basic,'tbriv')
        
        self.mode_radio_basic.show()
        self.mode_radio_tbr.show()
        self.mode_radio_tbriv.show()
        
        self.mode_box.pack_start(self.mode_radio_basic)
        self.mode_box.pack_start(self.mode_radio_tbr)
        self.mode_box.pack_start(self.mode_radio_tbriv)

        self.mode_box.show()
        
        self.big_box.pack_start(self.mode_box)

#makes the trace button

        self.trace_button = gtk.Button("Trace")
        self.trace_button.connect("clicked", self.trace, None)
        self.trace_button.show()
        self.big_box.pack_start(self.trace_button)

        #makes the graph buttons

        self.graph_box = gtk.HBox(False,10)

        self.graph_button_tbr = gtk.Button("Tbr Graph")
        self.graph_button_tbr.connect("clicked", self.tbr, None)
        self.graph_button_tbr.set_sensitive(False)
        self.graph_button_tbr.show()
  
        self.graph_box.pack_start(self.graph_button_tbr)

        self.graph_button_tbriv = gtk.Button("TbrIV Graph")
        self.graph_button_tbriv.connect("clicked", self.tbriv, None)
        self.graph_button_tbriv.set_sensitive(False)
        self.graph_button_tbriv.show()

        self.graph_box.pack_start(self.graph_button_tbriv)

        self.graph_button_traj = gtk.Button("Show Trajectories")
        self.graph_button_traj.connect("clicked", self.trajectory, None)
        self.graph_button_traj.set_sensitive(False)
        self.graph_button_traj.show()

        self.graph_box.pack_start(self.graph_button_traj)

        self.graph_box.show()

        self.big_box.pack_start(self.graph_box)

        

        #Make the sun

        self.sun_box = gtk.VBox(False,10)

        self.sun_drawing = gtk.DrawingArea()
        self.sun_drawing.set_size_request(400,400)
        

        self.sun_drawing.connect("expose-event", self.sunExpose)
        self.sun_drawing.connect("motion_notify_event", self.mouseMotion)
        self.sun_drawing.connect("button_press_event", self.mouseClicked)
        self.sun_drawing.connect("button_release_event", self.mouseReleased)
        self.sun_drawing.connect("drag_motion", self.mouseClicked)
        self.sun_drawing.connect("leave_notify_event", self.mouseExit)
        #self.sun_drawing.connect("button-release-event", self.sunExpose)

        self.sun_drawing.set_events(gtk.gdk.POINTER_MOTION_MASK |
                                    gtk.gdk.LEAVE_NOTIFY_MASK |
                                    gtk.gdk.BUTTON_PRESS_MASK |
                                    gtk.gdk.BUTTON_RELEASE_MASK)


        self.sun_drawing.show()

        self.sun_box.pack_start(self.sun_drawing)

    #makes the trjaectory buttons

        self.traj_box = gtk.HBox(False,10)

        self.traj_button_clear = gtk.Button("Clear")
        self.traj_button_clear.connect("clicked", self.clear, None)
        self.traj_button_clear.show()
  
        self.traj_box.pack_start(self.traj_button_clear)

        self.traj_button_select = gtk.Button("Select All")
        self.traj_button_select.connect('clicked', self.selectAll, None)
        self.traj_button_select.show()

        self.traj_box.pack_start(self.traj_button_select)

        self.streamer_button = gtk.Button("Add streamers")
        self.streamer_button.connect('clicked', self.streamers, None)
        self.streamer_button.show()

        self.traj_box.pack_start(self.streamer_button)

        self.traj_box.show()

        self.sun_box.pack_start(self.traj_box)

        self.sun_box.show()

#Add the menu bar

        self.menu_bar = gtk.MenuBar()

        self.menu_file = gtk.MenuItem("File")
        self.menu_file.show()

        self.menu_file_menu = gtk.Menu()
        self.menu_file_menu.show()

        self.menu_save = gtk.MenuItem("Save")
        self.menu_save.connect('activate',self.save)
        self.menu_save.show()

        self.menu_file_menu.append(self.menu_save)

        self.menu_batch = gtk.MenuItem("Make Batch")
        self.menu_batch.connect('activate', self.makeBatch)
        self.menu_batch.show()

        self.menu_file_menu.append(self.menu_batch)

        self.menu_exit = gtk.MenuItem("Exit")
        self.menu_exit.connect('activate', self.destroy)
        self.menu_exit.show()

        self.menu_file_menu.append(self.menu_exit)

        self.menu_file.set_submenu(self.menu_file_menu)

        self.menu_bar.append(self.menu_file)
        self.menu_bar.show()

        

        #Add the big_box to the window
        self.menu_box = gtk.VBox(False,1)
        
        self.big_big_box = gtk.HBox(False,10)
        self.big_big_box.pack_start(self.big_box)
        self.big_big_box.pack_start(self.sun_box)

        self.big_box.show()

        self.menu_box.add(self.menu_bar)
        self.menu_box.add(self.big_big_box)

        self.menu_box.show()
        
        self.window.add(self.menu_box)
        self.big_big_box.show()


    
        # and the window
        self.window.show()

        self.sun_drawing.window.resize(400,400)

        self.window.set_resizable(False)
        
    def main(self):
        # All PyGTK applications must have a gtk.main(). Control ends here
        # and waits for an event to occur (like a key press or mouse event).
        gtk.main()

# If the program is run directly or passed as an argument to the python
# interpreter then create a HelloWorld instance and show it


class StreamerWindow():

    def delete_event(self, widget, event, data=None):
        # If you return FALSE in the "delete_event" signal handler,
        # GTK will emit the "destroy" signal. Returning TRUE means
        # you don't want the window to be destroyed.
        # This is useful for popping up 'are you sure you want to quit?'
        # type dialogs.
        print "delete event occurred"

        # Change FALSE to TRUE and the main window will not be destroyed
        # with a "delete_event".
        return False

    def destroy(self, widget, data=None):
        print "destroy signal occurred"
        gtk.main_quit()

    def streamerExpose(self, widget, data=None):

        
        style = self.streamer_drawing.get_style()
        gc = style.text_aa_gc[gtk.STATE_NORMAL]

        gc.set_rgb_fg_color(gtk.gdk.Color('#ff0'))

        width = self.streamer_drawing.size_request()[0]
        height = self.streamer_drawing.size_request()[1]

        x = self.streamer_drawing.get_pointer()[0]
        y = self.streamer_drawing.get_pointer()[1]

        real_x = (x/float(width)-.5)*2
        real_y = (y/float(height)-.5)*2
        
        real_z = sqrt(1.0-real_x**2-real_y**2)

        if(math.isnan(real_z)):
            real_z = 0

        real_phi = math.acos(real_z/(sqrt(real_z**2+real_x**2)))
        real_theta = math.atan(real_y/(sqrt(real_z**2+real_x**2))) + math.pi/2.0

        if(real_x<0):
            real_phi = real_phi*-1
        

        self.layout = self.streamer_drawing.create_pango_layout('phi: '+str(real_phi)+' theta: '+str(real_theta))


        self.streamer_drawing.window.draw_arc(gc,True, 0,0,width-1,height-1,0,360*64)  

        gc.set_rgb_fg_color(gtk.gdk.Color('#000'))

        self.streamer_drawing.window.draw_arc(gc,False, 0,0,width-1,height-1,0,360*64)



        for stream in self.streamers:
            phi_sep = self.streamer_seperation*sin(stream[2])
            theta_sep = self.streamer_seperation*cos(stream[2])

            theta1 = stream[0] - theta_sep
            phi1 = stream[1] - phi_sep
            x1 = int((math.sin(theta1)*math.sin(phi1)+1)*width/2.0)
            y1 = int((math.cos(theta1)-1)*height/-2.0)

            theta2 = stream[0] + theta_sep
            phi2 = stream[1] + phi_sep
            x2 = int((math.sin(theta2)*math.sin(phi2)+1)*width/2.0)
            y2 = int((math.cos(theta2)-1)*height/-2.0)
            
            self.streamer_drawing.window.draw_arc(gc,True,
                                                  x1,
                                                  y1,
                                                  int(self.streamer_size*stream[3]),
                                                  int(self.streamer_size*stream[3]),0,360*64)
            self.streamer_drawing.window.draw_arc(gc,True,
                                                  x2,
                                                  y2,
                                                  int(self.streamer_size*stream[3]),
                                                  int(self.streamer_size*stream[3]),0,360*64)
            
        gc.set_rgb_fg_color(gtk.gdk.Color('#f00'))

        if(x>0 and y>0):
            phi_sep = self.streamer_seperation*sin(self.theta)
            theta_sep = self.streamer_seperation*cos(self.theta)

            theta1 = real_theta - theta_sep
            phi1 = real_phi - phi_sep
            x1 = int((math.sin(theta1)*math.sin(phi1)+1)*width/2.0)
            y1 = int((math.cos(theta1)-1)*height/-2.0)

            theta2 = real_theta + theta_sep
            phi2 = real_phi + phi_sep
            x2 = int((math.sin(theta2)*math.sin(phi2)+1)*width/2.0)
            y2 = int((math.cos(theta2)-1)*height/-2.0)
            
            self.streamer_drawing.window.draw_arc(gc,True,
                                                  x1,
                                                  y1,
                                                  int(self.streamer_size*real_z),
                                                  int(self.streamer_size*real_z),0,360*64)
            self.streamer_drawing.window.draw_arc(gc,True,
                                                  x2,
                                                  y2,
                                                  int(self.streamer_size*real_z),
                                                  int(self.streamer_size*real_z),0,360*64)

        if(self.released==1):
            self.streamers.add((real_theta,real_phi,self.theta,real_z))
            self.released = 0
            
        gc.set_rgb_fg_color(gtk.gdk.Color('#000'))

        self.streamer_drawing.window.draw_layout(gc, x=0, y=0, layout=self.layout)


        

    def mouseMotion(self, widget, data=None):
        self.streamer_drawing.queue_draw()

    def mouseExit(self, widget, data=None):
        self.streamer_drawing.queue_draw()

    def mouseScroll(self, widget, data=None):
        if(data.direction==gtk.gdk.SCROLL_DOWN):
            self.theta = self.theta + math.pi/12
        if(data.direction==gtk.gdk.SCROLL_UP):
            self.theta = self.theta - math.pi/12
        self.streamer_drawing.queue_draw()

    def mouseReleased(self, widget, data=None):
        self.released = 1
        self.streamer_drawing.queue_draw()

    def getStreamers(self):
        if(self.response == gtk.RESPONSE_ACCEPT):
            return self.streamers
        
    def __init__(self, parent, streamers):
        # create a new window

        if(streamers != None):
            self.streamers = streamers
        else:
            self.streamers = set()
        

        self.streamer_size = 16
        self.streamer_seperation = pi/12.0
        self.theta = 0
        self.released = 0

        self.dialog = gtk.Dialog("Streamers",
                                 parent,
                                 gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,
                                 (gtk.STOCK_CANCEL, gtk.RESPONSE_REJECT,
                                  gtk.STOCK_OK, gtk.RESPONSE_ACCEPT))

        self.streamer_drawing = gtk.DrawingArea()
        self.streamer_drawing.set_size_request(400,400)

        self.streamer_drawing.connect("expose-event", self.streamerExpose)

        self.streamer_drawing.connect("motion_notify_event", self.mouseMotion)
        #self.sun_drawing.connect("button_press_event", self.mouseClicked)
        self.streamer_drawing.connect("button_release_event", self.mouseReleased)
        #self.sun_drawing.connect("drag_motion", self.mouseClicked)
        self.streamer_drawing.connect("leave_notify_event", self.mouseExit)
        self.streamer_drawing.connect("scroll-event", self.mouseScroll)
        
        self.streamer_drawing.set_events(gtk.gdk.POINTER_MOTION_MASK |
                                         gtk.gdk.LEAVE_NOTIFY_MASK |
                                         gtk.gdk.BUTTON_PRESS_MASK |
                                         gtk.gdk.BUTTON_RELEASE_MASK |
                                         gtk.gdk.SCROLL_MASK)

        self.dialog.vbox.pack_start(self.streamer_drawing)
        self.streamer_drawing.show()

        self.dialog.set_resizable(False)

        self.response = self.dialog.run()

        if(self.response==gtk.RESPONSE_REJECT):
            self.dialog.destroy()

        if(self.response==gtk.RESPONSE_ACCEPT):
            self.dialog.destroy()

main = MainWindow()
main.main()

