import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdate
from parse_sol import parse_sol

COLOR1 = (77/256, 133/256, 189/256)
COLOR2 = (247/256, 144/256, 61/256)
COLOR3 = (89/256, 169/256, 90/256)


def dop_plot(path, data):
    t = []
    pdop = []
    total_sat_num = []
    use_sat_num = []

    for row in data:
        i = 0
        for l in row['epoch']:
            i = i + 1
            if i == 1:
                t.append(l)
            elif i == 6:
                total_sat_num.append(l)
            elif i == 7:
                use_sat_num.append(l)
            elif i == 8:
                pdop.append(l)
            else:
                continue

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))

    plt.plot_date(t, total_sat_num, '-', color=COLOR1, label='total')
    plt.plot_date(t, use_sat_num, '-', color=COLOR2, label='used')
    plt.legend(loc='best')

    plt.xlabel('Time [h]')
    plt.ylabel('Satellite Number')

    ax1 = ax.twinx()
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot_date(t, pdop, '-', color=COLOR3, label='PDOP')
    ax1.tick_params(axis='y', colors=COLOR3)
    plt.ylabel('PDOP')

    plt.show()
    save_name = path[-8:-4] + '_dop.png'
    save_dir = path[:-9]
    save_path = os.path.join(save_dir, save_name)
    fig.savefig(save_path, bbox_inches='tight', dpi=300)


def enu_plot(path, data):
    t = []
    e = []
    n = []
    u = []
    for row in data:
        i = 0
        for l in row['epoch']:
            if row['epoch'][4] == 'NONE':
                break
            i = i + 1
            if i == 1:
                t.append(l)
            else:
                continue

    for row in data:
        i = 0
        for l in row['pos']:
            i = i + 1
            if i == 1:
                e.append(l)
            elif i == 2:
                n.append(l)
            elif i == 3:
                u.append(l)
            else:
                continue

    if e.__len__() <= 0 or n.__len__() <= 0 or u.__len__() <= 0 or t.__len__() <= 0:
        return

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    plt.plot_date(t, e, color=COLOR1, ms=0.5, label='east')
    plt.plot_date(t, n, color=COLOR2, ms=0.5, label='north')
    plt.plot_date(t, u, color=COLOR3, ms=0.5, label='up')
    ax.axhline(y=0.1, color='y', linestyle='-', linewidth=0.4)
    ax.axhline(y=-0.1, color='y', linestyle='-', linewidth=0.4)
    plt.ylim([-0.5, 0.5])
    plt.legend()
    plt.xlabel('Time [h]')
    plt.ylabel('Position Error [m]')
    # plt.show()

    save_name = path[-8:-4] + '_pos.png'
    save_dir = path[:-9]
    save_path = os.path.join(save_dir, save_name)
    fig.savefig(save_path, bbox_inches='tight', dpi=300)


def plot_isb():
    a = 1


def trp_plot(path, data):
    t = []
    dry = []
    wet = []
    for row in data:
        i = 0
        for l in row['epoch']:
            i = i + 1
            if i == 1:
                t.append(l)
            else:
                continue

    for row in data:
        i = 0
        for l in row['trp']:
            i = i + 1
            if i == 1:
                dry.append(l)
            elif i == 2:
                wet.append(l)
            else:
                continue

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    plt.plot_date(t, dry, '-', color=COLOR1)

    plt.xlabel('Time [h]')
    plt.ylabel('Dry [m]')

    ax1 = ax.twinx()
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot_date(t, wet, '-', color=COLOR3)
    ax1.tick_params(axis='y', colors=COLOR3)
    plt.ylabel('Wet [m]')

    plt.show()
    save_name = path[-8:-4] + '_trp.png'
    save_dir = path[:-9]
    save_path = os.path.join(save_dir, save_name)
    fig.savefig(save_path, bbox_inches='tight', dpi=300)


def ion_plot(path, data, sat_id):
    t = []
    ion = []
    for row in data:
        i = 0
        for l in row['epoch']:
            i = i + 1
            if i == 1:
                t.append(l)
            else:
                continue

    for row in data:
        for l in row['sat']:
            if sat_id == l[0]:
                j = 0
                for s in l:
                    j = j + 1
                    if j == 19:
                        ion.append(s)
                    else:
                        continue
            else:
                continue

    ions = []
    ts = []
    for i in range(0, ion.__len__()):
        if ion[i] is not None:
            ions.append(ion[i])
            ts.append(t[i])

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    # ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    # ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    # ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    plt.plot_date(ts, ions)
    plt.show()


DOP_PLOT = 1
ENU_PLOT = 1
ISB_PLOT = 0
TRP_PLOT = 0
ION_PLOT = 0
AMB_PLOT = 0
SAT_ID = 'G10'
if __name__ == '__main__':
    sol_path = '../../example/result_wum/PPP_STATIC/G_DF_IF/faa1.pos'
    data = parse_sol(sol_path)
    if ENU_PLOT:
        enu_plot(sol_path, data)
    if DOP_PLOT:
        dop_plot(sol_path, data)
    if ISB_PLOT:
        a = 1
    if TRP_PLOT:
        trp_plot(sol_path, data)
    if ION_PLOT:
        ion_plot(sol_path, data, SAT_ID)

