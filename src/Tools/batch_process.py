import subprocess
import datetime
from dateutil.relativedelta import relativedelta
import platform


PRC_MODE = ["PPP-STATIC"]
PRC_SYS = ["G", "GB"]
PRC_FRQ = [2]
PRC_ION = [3, 4]  # 3 if 4 uc


def ppplib(argv):
    if platform.system() == 'Windows':
        cmd = 'C:\\Users\\chenc\\Desktop\\PPPLib_v1.0\\bin/PPPMain '+argv
    else:
        cmd = '/home/cc/Codings/PPPLib_v1.0/bin/PPPMain '+argv
    print(cmd)
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    time_start = datetime.datetime(2019, 12, 8, 0, 0, 0)
    time_end = datetime.datetime(2019, 12, 8, 0, 0, 0)
    time_idx = time_start

    while 1:
        if (time_idx - time_end).total_seconds() > 0:
            break

        prc_opt = "-do 1"
        prc_date = "-pd " + str(time_idx.year) + '/' + str(time_idx.month) + '/' + str(time_idx.day)
        for i in PRC_MODE:
            prc_mod = "-md " + i
            for j in PRC_SYS:
                prc_sys = "-sys " + j
                for k in PRC_ION:
                    prc_ion = "-ion " + str(k)
                    for n in PRC_FRQ:
                        prc_frq = "-frq " + str(n)
                        argv = prc_opt + " " + prc_date + " " + prc_mod + " " + prc_sys + " " + prc_ion + " " + prc_frq + " "+ \
                               "-level 128"
                        ppplib(argv)
                        print(argv)

        # dir = "/home/cc/test_data/2019/208/335/result_wum/PPP_KINE"
        # ps.walk_files(dir, [], " ")

        time_idx = time_idx + relativedelta(days=1)