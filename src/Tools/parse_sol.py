# coding: utf-8
from collections import namedtuple
from datetime import datetime


class Site:
    def __init__(self):
        self.station = ''
        self.rec_type = ''
        self.ant_type = ''
        self.sta_pos = []
        self.nav_sys = ''
        self.trp_mode = ''
        self.ion_mode = ''
        self.sat_eph = ''
        self.epoch = []
        self.pos = ''
        self.clk = ''
        self.trp = []
        self.sat = []


site = Site()
data = []


def parse_sol(path):
    with open(path, 'r', encoding='utf-8') as f:
        n = 0
        t = 0
        for line in f:
            n += 1
            line = line.strip()
            if n < 16:
                if n == 1:
                    if line.startswith('+'):
                        if line == '+ PPPLib HEADER':
                            pass
                        else:
                            raise ValueError
                    else:
                        raise ValueError

                elif n == 2:
                    if line.startswith('Station :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.station = __.strip()

                elif n == 3:
                    if line.startswith('Rec Type:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.rec_type = __.strip()

                elif n == 4:
                    if line.startswith('Ant Type:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.ant_type = __.strip()

                elif n == 5:
                    if line.startswith('Sta  Pos:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    __ = [float(i) for i in __.split(' ') if i.strip()]
                    site.sta_pos = __

                elif n == 6:
                    if line.startswith('Nav  Sys:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.nav_sys = __.strip()

                elif n == 7:
                    if line.startswith('Trp Mode:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.trp_mode = __.strip()

                elif n == 8:
                    if line.startswith('Ion Mode:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.ion_mode = __.strip()

                elif n == 9:
                    if line.startswith('Sat  Eph:'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.sat_eph = __.strip()

                elif n == 10:
                    if line.startswith('$EPOCH  :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    __ = [i.strip() for i in __.split('|')]
                    site.epoch = __

                elif n == 11:
                    if line.startswith('$POS    :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.pos = __.strip()

                elif n == 12:
                    if line.startswith('$CLK    :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.clk = __.strip()

                elif n == 13:
                    if line.startswith('$TRP    :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.trp = [i.strip() for i in __.split('|')]

                elif n == 14:
                    if line.startswith('$SAT    :'):
                        pass
                    else:
                        raise ValueError
                    _, __ = line.split(':', 1)
                    site.sat = [i.strip() for i in __.split('|')]

                if n == 15:
                    if line.startswith('-'):
                        if line == '- PPPLib HEADER':
                            pass
                        else:
                            raise ValueError
                    else:
                        raise ValueError

                continue
            # 空行 新建容器，新建数据容器
            if not line:
                data.append({
                    'epoch': [],
                    'pos': [],
                    'clk': [],
                    'trp': [],
                    'sat': [],
                })
                if t > 0:
                    print(n, t, len(data[t-1]['sat']))
                    # if len(data[t-1]['sat']) != 10:
                    #     raise ValueError
                t += 1
            if line .startswith('$EPOCH'):
                tmp = data[t-1]

                line = [i.strip() for i in line.split(' ') if i.strip()]
                if len(line) != 11:
                    raise ValueError
                time = ' '.join(line[1:3])
                time = datetime.strptime(time, '%Y/%m/%d %H:%M:%S.0')
                line[3] = int(line[3])  # WEEK
                line[4] = float(line[4])  # WOS
                line[5] = int(line[5])  # Epoch
                line[7] = int(line[7])  # ObsNum
                line[8] = int(line[8])  # ValidSat
                line[9] = float(line[9])  # PDOP
                line[10] = float(line[10])  # SIGMA0

                tmp['epoch'] = [time] + line[3:]
            elif line.startswith('$POS'):
                tmp = data[t - 1]
                line = [float(i.strip()) for i in line.split(' ')[1:] if i.strip()]
                if len(line) != 3:
                    raise ValueError
                tmp['pos'] = line

            elif line.startswith('$CLK'):
                tmp = data[t - 1]
                line = [float(i.strip()) for i in line.split(' ')[1:] if i.strip()]
                if len(line) != 6:
                    raise ValueError
                tmp['clk'] = line

            elif line.startswith('$TRP'):
                tmp = data[t - 1]
                line = [float(i.strip()) for i in line.split(' ')[1:] if i.strip()]
                if len(line) != 2:
                    raise ValueError
                tmp['trp'] = line

            elif line.startswith('$SAT'):
                tmp = data[t - 1]
                line = [i.strip() for i in line.split(' ')[1:] if i.strip()]
                if len(line) != 22:
                    raise ValueError
                tmp['sat'].append(line)

    # 获取不重复的 sat id
    sat_id_set = set()
    for row in data:
        for l in row['sat']:
            sat_id_set.add(l[0])

    print('获取不重复的 sat id', len(sat_id_set))

    # 补齐每个历元缺失的 sat
    for row in data:
        _ls = [l[0] for l in row['sat']]
        for i in sat_id_set - set(_ls):
            row['sat'].append([i] + [None] * 21)
        row['sat'].sort(key=lambda x: x[0])

    for row in data:
        print(len(row['sat']))

    return data


fp = '../../example/result_wum/PPP_KINE/G_DF_IF/jfng.pos'
if __name__ == '__main__':
    data = parse_sol(fp)



