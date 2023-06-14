import numpy as np
from bulk_generator import *


class Taper_wg:
    def __init__(self, air_width=1000, length_limit=2000000, quantization=50, narrow_end_position=(0, 0),
                 narrow_end_width=400) -> None:
        self.air_width = air_width
        self.length_limit = length_limit
        self.quantization = quantization
        self.narrow_end_position = narrow_end_position
        self.narrow_end_width = narrow_end_width

    def parallelogram_from_left_bottom(self, l, w, dy):
        return np.array(
            [[0, 0],
             [l, dy],
             [l, dy + w],
             [0, w]
             ]
        )

    def linear(self, start, dx, dy, air_width=None):
        aw = air_width if air_width is not None else self.air_width
        upper = shift(self.parallelogram_from_left_bottom(dx, aw, dy), start)
        end = shift(start, (dx, dy))
        return upper, end

    def sin(self, start, dx, dy, air_width=None):
        aw = air_width if air_width is not None else self.air_width
        k = np.pi / dx / 2
        num_segments = int(dx / self.quantization)  # Number of linear segments to approximate the sine curve
        segment_length = dx / num_segments
        A = (-dy) / (np.sin(np.pi / 2 - k * dx) - 1)

        prev_end = start.tolist()
        structures = [prev_end]
        for i in range(num_segments):
            x_start = prev_end[0]
            x_end = x_start + segment_length
            y_end = A * np.sin(k * (x_end - start[0] - dx) + np.pi / 2) + dy - A + start[1]
            end = [x_end, y_end]
            structures.append(end)
            prev_end = end
        original_structures = structures.copy()
        structures.reverse()
        original_structures += [[var[0], var[1] + aw] for var in structures]

        return np.array(original_structures), np.array(prev_end)

    def make_many(self, *args):
        structures = []
        start = self.narrow_end_position

        if not isinstance(start, np.ndarray):
            start = np.array(start)

        for lst in args:
            if len(lst) == 4:
                the_string, dx, dy, aw = lst
            else:
                the_string, dx, dy = lst
                aw = self.air_width
            if the_string in ["linear", "l", "lin"]:
                the_structure, end = self.linear(start, dx, dy, aw)
            elif the_string in ["sine", "sin", "s"]:
                the_structure, end = self.sin(start, dx, dy, aw)
            structures.append(the_structure)
            start = end
        structures = shift(structures, (0, self.narrow_end_width / 2))
        structures += flip_points(structures, [(0, 0), (1, 0)])
        return structures

    def preview(self, var):
        plot_many(lst_polyline=var, show_coords=True)


# test=Taper_wg()
# a=test.make_many(["l",4000,0],["s",10000,1000],["l",1000,0],["s",15000,1500],["l",1000,0],["l",40000,1500],["l",100000,0])


connect_wire_w = 1000
connect_wire_l = 12500

taper1_butt = 2000
taper1_shrink = 1250 * 2
taper1_l = 15500
taper1_linear_part_l = 1200
N = 20

tapper_connector_l = 4520

taper2_butt = 2000
taper2_shrink = 1100 * 2
taper2_l = 15500
taper2_linear_part_l = 1200

taper3_l = 40000
taper3_shrink = 1400 * 2

siwg_w = 8000
siwg_l = 200000
siwg_n = 12


def connect_wire():
    return np.array([[0, 0],
                     [connect_wire_w, 0],
                     [connect_wire_w, connect_wire_l],
                     [0, connect_wire_l]]
                    )


# taper3_shrink=(siwg_w-taper2_shrink-taper1_shrink)
def taper1():
    lst = []
    lst.append([0, 0])
    lst.append([0, taper1_linear_part_l])
    for i in range(N):
        y = taper1_linear_part_l + (taper1_l - taper1_linear_part_l) * ((i + 1) / N)
        x = taper1_shrink / 2 * np.sin(0.5 * np.pi * ((i + 1) / N))
        lst.append([x, y])
    lst.append([x + taper1_butt, y])
    lst.append([x + taper1_butt, 0])
    return np.array(lst)


def taper2():
    lst = []
    lst.append([0, 0])
    lst.append([0, taper2_linear_part_l])
    for i in range(N):
        y = taper2_linear_part_l + (taper2_l - taper2_linear_part_l) * ((i + 1) / N)
        x = taper2_shrink / 2 * np.sin(0.5 * np.pi * ((i + 1) / N))
        lst.append([x, y])
    lst.append([x + taper2_butt, y])
    lst.append([x + taper2_butt, 0])
    return np.array(lst)


def taper_connector():
    return np.array([[0, 0],
                     [taper1_butt, 0],
                     [taper1_butt, tapper_connector_l],
                     [0, tapper_connector_l]]
                    )


def taper3():
    return np.array([[0, 0],
                     [taper2_butt, 0],
                     [taper2_butt + taper3_shrink / 2, taper3_l],
                     [taper3_shrink / 2, taper3_l]]
                    )


def SiWG():
    result = []
    standard_shape = np.array([[0, 0],
                               [taper2_butt, 0],
                               [taper2_butt, siwg_l],
                               [0, siwg_l]]
                              )
    standard_shape = shift(standard_shape, (taper1_shrink / 2 + taper2_shrink / 2 + taper3_shrink / 2,
                                            connect_wire_l + taper1_l + tapper_connector_l + taper2_l + taper3_l))
    for i in range(siwg_n):
        result.append(standard_shape)
        standard_shape = shift(standard_shape, (0, siwg_l))
    return result


def Yoda_SiWG_right():
    wire = connect_wire()
    taper_1 = taper1()
    taper_1 = shift(taper_1, (0, connect_wire_l))
    t_connector = taper_connector()
    t_connector = shift(t_connector, (taper1_shrink / 2, connect_wire_l + taper1_l))
    taper_2 = taper2()
    taper_2 = shift(taper_2, (taper1_shrink / 2, connect_wire_l + taper1_l + tapper_connector_l))
    taper_3 = taper3()
    taper_3 = shift(taper_3,
                    (taper1_shrink / 2 + taper2_shrink / 2, connect_wire_l + taper1_l + tapper_connector_l + taper2_l))
    straight_part = SiWG()
    result = [wire, taper_1, t_connector, taper_2, taper_3] + straight_part
    # plot_many([],result)
    return result


def Yoda_SiWG_left():
    r = Yoda_SiWG_right()
    l = flip_points(r, [(0, 0), (0, 1)])
    return r
