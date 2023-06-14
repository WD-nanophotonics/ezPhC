import ezdxf
import numpy as np
from typing import List, Tuple
from typing import Literal
from bulk_generator import *

from SiWG import Yoda_SiWG_right,Yoda_SiWG_left


class Lattice:  # the main object
    def __init__(self, period: float, outline: List[Tuple[float, float]], orientation: float,
                 lattice_type: str) -> None:
        self.period = period  # lattice constant
        self.outline = [(self.period*x,self.period*y) for x,y in outline]  # an outline to define the range of the PhC
        self.orientation = orientation  # lattice direction, GK and GX by default
        self.lattice_type = lattice_type  # supports triangular, square and honeycomb
        self.points = self.get_points()  # the lattice positions. A list storing some (N,2)-shaped numpy arrays

    def FreePattern(self, *args: List[Tuple[float, float]]):  # put a customized pattern at each position in the lattice
        return self.inner1_freepattern(self, *args)

    def PolygonPattern(self, n, r, theta, n2=None, r2=None, theta2=None): # use this convenient predifined method if your pattern is polygon!!!
        pattern1 = [
            (r * np.sin(2 * i * np.pi / n + np.deg2rad(theta)), r * np.cos(2 * i * np.pi / n + np.deg2rad(theta))) for i
            in range(n)]
        if (n2 != None) * (r2 != None) * (theta2 != None):  # for honeycomb lattice, define patterns for a second (sub)lattice
            pattern2 = [(r2 * np.sin(2 * i * np.pi / n2 + np.deg2rad(theta2)),
                         r2 * np.cos(2 * i * np.pi / n2 + np.deg2rad(theta2))) for i in range(n2)]
            return self.inner1_freepattern(self, pattern1, pattern2)
        else:
            return self.inner1_freepattern(self, pattern1)

    def get_points(self): # the main function of Lattice class. Builds the actual lattice configuration
        poly = np.array(self.outline)
        shift_vector = np.mean(poly, axis=0)
        shifted_poly = poly - shift_vector
        large_enough_N = 5 * int(np.amax(shifted_poly)/self.period)  #create a large enough bulk lattice and later cut it off
        if self.lattice_type in ["triangular", "tri", "t"]:
            pos = maketriangularlattice(period=self.period, size=large_enough_N, draw=False, show_coords=False,
                                        angle=self.orientation, shift=shift_vector)
            pos = cut(pos, poly)
            return [pos]
        elif self.lattice_type in ["square", "sqr", "s"]:
            pos = makesquarelattice(period=self.period, size=large_enough_N, draw=False, show_coords=False,
                                    angle=self.orientation, shift=shift_vector)
            pos = cut(pos, poly)
            return [pos]
        elif self.lattice_type in ["honeycomb", "hon", "hc", "h"]: # HC lattice has 2 sublattices, requires separate operation
            pos1 = maketriangularlattice(period=self.period, size=large_enough_N, draw=False, show_coords=False)
            pos2 = maketriangularlattice(period=self.period, size=large_enough_N, draw=False, show_coords=False,
                                         shift=(self.period / 2, self.period / 2 / np.sqrt(3)))

            def separate_operation(pos):
                pos = rotate(pos, self.orientation)
                pos = shift(pos, shift_vector)
                pos = cut(pos, poly)
                return pos

            pos1 = separate_operation(pos1)
            pos2 = separate_operation(pos2)
            return [pos1, pos2]


    def preview_lattice(self, show_outline=False): #call pyplot
        pos = self.points

        if show_outline:
            plot_many(pos, [self.outline], show_coords=True)
        else:
            plot_many(pos, [], show_coords=True)

    def align_to_outline(self, main_index: Literal["top", "bottom", "left", "right"],  # align the lattice to the outline corners
                         secondary_index: Literal["top", "bottom", "left", "right"], sublattice=None):
        pos = self.points[0]
        if self.lattice_type in ["honeycomb", "hon", "hc", "h"] and sublattice == 2:
            pos = self.points[1]

        point_to_move = find_point(pos, main_index, secondary_index)
        target_position = find_point(np.array(self.outline), main_index, secondary_index)
        shift_vector = target_position - point_to_move

        self.points = [shift(var, shift_vector) for var in self.points]

    def shift(self, vector):
        self.points=shift(self.points,vector)

    def rotate(self,angle,base_point):
        self.points=rotate(self.points,angle,base_point)

    def mirror(self, axis, keep_original=True, glide=False):

        reflected =flip_points(self.points, axis)

        # print(reflected)
        if glide:
            reflected = shift_along_axis(reflected, axis, self.period / 2)

        if keep_original:

            self.points = [np.concatenate((self.points[i], reflected[i]), axis=0) for i in range(len(self.points))]

        else:
            self.points = reflected

    class inner1_freepattern:
        def __init__(self, outer_instance, *args: List[Tuple[float, float]]) -> None:
            self.outer_instance = outer_instance
            self.args = [var for var in args]
            self.pattern = self.get_pattern()

        def get_pattern(self):
            lattice_grids = self.outer_instance.points
            return [[self.args[i] + P for P in lattice_grids[i]] for i in range(len(self.args))]

        """def mirror_(self, axis, keep_original=True, glide=False):
            result = []

            for sublattice in self.pattern:
                sub = []
                for pats in sublattice:
                    pats = flip_points(pats, axis=axis)
                    if glide:
                        pats = shift_along_axis(pats, axis, self.outer_instance.period / 2)
                    sub.append(pats)
                result.append(sub)
                if keep_original:
                    result.append(sublattice)
            self.pattern = result"""

        def mirror(self, axis, keep_original=True, glide=False):
            reflected = flip_points(self.pattern, axis)
            # print(reflected)
            if glide:
                reflected = shift_along_axis(reflected, axis, self.outer_instance.period / 2)
            if keep_original:

                self.pattern = [self.pattern[i]+reflected[i] for i in range(len(self.pattern))]

            else:
                self.pattern = reflected

        def truncate(self,polyline): # clumsy style, needs revision
            result = []
            polyline=[(self.outer_instance.period*x,self.outer_instance.period*y) for x,y in polyline]
            for sublattice in self.pattern:
                sub = []
                for pats in sublattice:
                    #print("patterns=",pats)
                    center=np.mean(pats,axis=0)
                    #print("c=",center)
                    if point_in_polygon(tuple(center),polyline):
                        sub.append(pats)
                result.append(sub)
            self.pattern = result

        def add_pattern(self,other_pattern):
            other_pattern=convet_to_two_layer_pattern_list(other_pattern)
            self.pattern+=other_pattern
        def preview_pattern(self, show_outline=False):
            pos = self.pattern

            pos = [inner_list for outer_list in pos for inner_list in outer_list]

            if show_outline:
                plot_many([], pos + [self.outer_instance.outline], show_coords=False)
            else:
                plot_many([], pos, show_coords=False)

        def save_to_dxf(self, fname, version, text=None, textposition=None,textsize=1.0):
            pos = self.pattern
            pos = [inner_list for outer_list in pos for inner_list in outer_list]
            add_polygons_to_dxf(pos, fname)
            if text != None:
                add_text_to_dxf(fname, text, textposition,textsize)
