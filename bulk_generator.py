import ezdxf
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import orient
from shapely.geometry.polygon import Polygon
from data_handler import *

def show_coorinate(fig,ax,scatter):
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = f"({pos[0]:.2f}, {pos[1]:.2f})"
        annot.set_text(text)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()


def plot_lattice(pos,show_coords=False):
    fig, ax = plt.subplots()
    scatter = ax.scatter(pos[:, 0], pos[:, 1])
    if show_coords:
        show_coorinate(fig,ax,scatter)
    plt.axis("equal")
    plt.show()


def plot_polygon(poly):
    # Extract x and y coordinates from each point in the list
    x, y = zip(*poly)

    # Close the polygon by adding the first point to the end
    x += (x[0],)
    y += (y[0],)

    # Plot the polygon using Matplotlib
    plt.plot(x, y, 'b-')
    plt.fill(x, y, alpha=0.2)
    plt.axis('equal')
    plt.show()


def point_in_polygon(pos,polygon):
    poly = Polygon(polygon)
    poly=orient(poly,sign=1)
    if poly.contains(Point(pos)):
        return True
    else:
        return False

def cut(points, polygon):
    """Remove points outside a polygon"""
    poly = Polygon(polygon)
    poly=orient(poly,sign=1)
    mask = np.array([Point(point).within(poly) for point in points])
    return points[mask]

def plot_many(lst_points=None, lst_polyline=None, show_coords=True):
    fig, ax = plt.subplots()
    scatter_holder = []
    annot_holder = []

    if lst_points is not None:
        for pos in lst_points:
            scatter = ax.scatter(pos[:, 0], pos[:, 1])
            scatter_holder.append(scatter)
            annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                 bbox=dict(boxstyle="round", fc="w"),
                                 arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            annot_holder.append(annot)

    if lst_polyline is not None:
        for poly in lst_polyline:
            x, y = zip(*poly)
            # Close the polygon by adding the first point to the end
            x += (x[0],)
            y += (y[0],)
            ax.plot(x, y, 'b-')

    def update_annot(ind, scatter):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot = annot_holder[scatter_holder.index(scatter)]
        annot.xy = pos
        text = f"({pos[0]:.2f}, {pos[1]:.2f})"
        annot.set_text(text)

    def hover(event):
        for scatter, annot in zip(scatter_holder, annot_holder):
            vis = annot.get_visible()
            if event.inaxes == scatter.axes:
                cont, ind = scatter.contains(event)
                if cont:
                    update_annot(ind, scatter)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

    if show_coords:
        fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.axis("equal")
    plt.show()



"""def plot_many2(lst_points=None,lst_polyline=None,show_coords=True):
    fig, ax = plt.subplots()
    if lst_points is not None:
        for pos in lst_points:
            scatter = ax.scatter(pos[:, 0], pos[:, 1])
    if lst_polyline is not None:
        for poly in lst_polyline:
            x, y = zip(*poly)
            # Close the polygon by adding the first point to the end
            x += (x[0],)
            y += (y[0],)
            ax.plot(x, y, 'b-')

    if show_coords:
        annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)

        def update_annot(ind):
            pos = scatter.get_offsets()[ind["ind"][0]]
            annot.xy = pos
            text = f"({pos[0]:.2f}, {pos[1]:.2f})"
            annot.set_text(text)

        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = scatter.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.axis("equal")
    plt.show()
"""
def shift(pos,vector): #takes numpy array as input
    if isinstance(pos,np.ndarray):
        return pos+vector

    else:
        result=[]
        for var in pos:
            new_var=shift(var,vector)
            result.append(new_var)
        return result

def rotate(pos,angle,origin=(0,0)):

    if isinstance(pos,np.ndarray):
        pos -= origin
        pos = np.dot(pos,
                     np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]]))  # rotate by angle
        pos += origin
        return pos
    else:
        result=[]
        for var in pos:
            new_var=rotate(var,angle)
            result.append(new_var)
        return result

def maketriangularlattice(period, size, draw=False, show_coords=True,angle=0.0,shift=(0,0)):
    rows = size
    cols = size
    pos = np.array([[(i+j%2/2)*period, j*np.sqrt(3)/2*period] for i in range(cols) for j in range(rows)])
    pos -= np.mean(pos, axis=0)  # center at the origin
    pos = np.dot(pos, np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]))  # rotate by angle
    # Shift the lattice by the given vector
    pos += np.array(shift)
    if draw:
        plot_lattice(pos,show_coords)

    return pos


def makesquarelattice(period, size, draw=False, show_coords=True,angle=0.0,shift=(0,0)):
    rows = size
    cols = size
    pos = np.array([[i*period, j*period] for i in range(cols) for j in range(rows)])
    pos=pos.astype(np.float64)
    pos -= np.mean(pos, axis=0)  # center at the origin
    pos = np.dot(pos, np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]))  # rotate by angle
    # Shift the lattice by the given vector
    pos += np.array(shift)
    if draw:
        plot_lattice(pos,show_coords)

    return pos

def find_point(points, main_index, secondary_index):  #find points in a numpy array in main-secondary order

    # Check that main_index and secondary_index are valid
    valid_indices = ["top", "bottom", "left", "right"]
    if main_index not in valid_indices or secondary_index not in valid_indices:
        raise ValueError("main_index and secondary_index must be elements of ['top', 'bottom', 'left', 'right']")

    # Determine which coordinate to use for the primary and secondary condition
    if main_index == "top":
        primary_coord = 1
        sort_order = -1
    elif main_index == "bottom":
        primary_coord = 1
        sort_order = 1
    elif main_index == "left":
        primary_coord = 0
        sort_order = 1
    elif main_index == "right":
        primary_coord = 0
        sort_order = -1

    # Find the point with the largest/smallest value of the primary coordinate

    sorted_points = points[np.argsort(points[:, primary_coord])[::sort_order]]
    primary_point = sorted_points[0]

    # Find the only point that satisfies the secondary condition
    if secondary_index == "top":
        secondary_coord = 1
        secondary_sort_order = -1
    elif secondary_index == "bottom":
        secondary_coord = 1
        secondary_sort_order = 1
    elif secondary_index == "left":
        secondary_coord = 0
        secondary_sort_order = 1
    elif secondary_index == "right":
        secondary_coord = 0
        secondary_sort_order = -1
    secondary_points = sorted_points[sorted_points[:, primary_coord] == primary_point[primary_coord]]
    sorted_points = secondary_points[np.argsort(secondary_points[:, secondary_coord])[::secondary_sort_order]]

    return sorted_points[0]


def find_largest_differences(array):
    x_values = array[:, 0]
    y_values = array[:, 1]
    max_x_diff = np.max(x_values) - np.min(x_values)
    max_y_diff = np.max(y_values) - np.min(y_values)
    return max_x_diff, max_y_diff

def flip_points(points, axis):
    """
    Flips a set of 2D points about a line defined by two points.
    """
    if isinstance(points,np.ndarray):
        # Convert the axis to a vector
        axis_vec = np.array(axis[1]) - np.array(axis[0])

        # Compute the normal to the axis
        normal_vec = np.array([-axis_vec[1], axis_vec[0]])

        # Compute the distance from the origin to the axis
        d = np.dot(normal_vec, np.array(axis[0]))

        # Compute the projection of the points onto the axis
        #proj = np.dot(points, axis_vec) / np.dot(axis_vec, axis_vec)

        # Compute the projection of the points onto the normal
        norm = np.dot(points, normal_vec) / np.dot(normal_vec, normal_vec)

        # Compute the reflected points
        reflected = points - 2 * ( norm[:, np.newaxis] * normal_vec)

        return reflected

    else:
        result=[]
        for var in points:
            new_var=flip_points(var,axis)
            result.append(new_var)
        return result


def shift_along_axis(point,axis,dist):
    if isinstance(point,np.ndarray):
        distance=sum([(axis[0][i]-axis[1][i])**2 for i in range(len(axis[0]))])**0.5
        v=[dist*(axis[1][i]-axis[0][i])/distance for i in range(len(axis[0]))]

        return point + v
    else:
        result = []
        for var in point:
            new_var = shift_along_axis(var, axis,dist)
            result.append(new_var)
        return result


"""def save_polygons_to_dxf(polygons, filename):

    doc = ezdxf.new(dxfversion='R2010')
    msp = doc.modelspace()
    if isinstance(polygons,np.ndarray):
        polyline = msp.add_polyline2d(polygons)
        polyline.close(True)
    else:
        for polygon in polygons:
            save_polygons_to_dxf(polygon,filename)

    doc.saveas(filename)

def save_polygons_to_dxf2(polygons, filename):

    doc = ezdxf.new(dxfversion='R2010')
    msp = doc.modelspace()
    if isinstance(polygons,np.ndarray):
        polyline = msp.add_polyline2d(polygons)
        polyline.close(True)
    else:
        for polygon in polygons:
            save_polygons_to_dxf(polygon,filename)

    doc.saveas(filename)"""

def add_polygons_to_dxf(polygons, filename):
    try:
        doc = ezdxf.readfile(filename)
    except:
        doc = ezdxf.new(dxfversion='R2010')
    msp = doc.modelspace()
    def write(polygons):

        if isinstance(polygons,np.ndarray):
            polyline = msp.add_polyline2d(polygons)
            polyline.close(True)

        else:
            for polygon in polygons:
                write(polygon)


    write(polygons)
    doc.saveas(filename)




def add_text_to_dxf(filename, text, position,size):
    # Load the DXF file
    doc = ezdxf.readfile(filename)
    msp = doc.modelspace()

    # Add the text
    msp.add_text(text, dxfattribs={'insert': position,"height":size})

    # Save the modified DXF file
    doc.saveas(filename)

"""points = np.array([[0, 3], [1, 1], [2, 4], [5, 2], [3, 0]])
print(find_point(points, "top", "left"))  # Output: [2, 4]"""
