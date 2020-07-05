import os
import platform
import ftplib
import datetime
from dateutil.relativedelta import relativedelta
import gzip
import unlzw
import subprocess
import cmn_func as common


LOCAL_DIR = 'E:\\Codings\\data_tmp\\'


class FtpDown:
    ftp = ftplib.FTP()

    def __init__(self, host, port=21):
        self.ftp.connect(host, port)

    def login(self):
        self.ftp.login()
        print(self.ftp.welcome)

    def file_classification(self, remote_dir, condition, file_type):
        """
        :param remote_dir:
        :param condition: sta_list or country_list for obs, ac_list for pre
        :param file_type: 'obs' or 'pre'
        :return:
        """
        files = []
        files_tmp = self.ftp.nlst(remote_dir)

        if file_type == 'obs':
            obs_files = [file for file in files_tmp if file.endswith('.crx.gz')]
            if condition:
                files = [file for file in obs_files if file[0:4] in condition or file[6:9] in condition]
            else:
                files = obs_files

        elif file_type == 'pre':
            if condition is 'GBM':
                pre_files = [file for file in files_tmp if (file.endswith('.CLK.gz') or file.endswith('.SP3.gz')
                                                            or file.endswith('ERP.gz')) and 'RAP' in file]
                for i in range(0, pre_files.__len__()-1):
                    files.append(pre_files[i][24:])
            else:
                pre_files = [file for file in files_tmp if (file.endswith('.CLK.gz') or file.endswith('.SP3.gz')
                                                            or file.endswith('ERP.gz')) and 'FIN' in file]
                if condition:
                    files = [file for file in pre_files if file[0:3] in condition]
        elif file_type == 'ion_dcb':
            files = [file for file in files_tmp if file.split('/')[-1] in condition]
        elif file_type == 'igs':
            files = [file for file in files_tmp if file in condition]

        return files

    def down_file_single(self, local_path, remote_file_name):
        if os.path.exists(local_path):
            print('uncompress file exist %s, please delete, re-download!' % local_path)
            return False
            # os.remove(local_path)
        file_handler = open(local_path, "wb")
        print(file_handler)

        self.ftp.retrbinary('RETR ' + remote_file_name, file_handler.write)
        file_handler.close()
        return True

    def down_file_batch(self, local_dir, remote_dir):
        print("RemoteDir:", remote_dir)
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)

        self.ftp.cwd(remote_dir)

        remote_files = self.ftp.nlst()
        for file in remote_files:
            local_path = os.path.join(local_dir, file)
            print(self.ftp.nlst(file))
            if file.find(".") == -1:
                if not os.path.exists(local_path):
                    os.makedirs(local_path)
                    self.down_file_single(local_path, file)
                else:
                    self.down_file_single(local_path, file)
        self.ftp.cwd("..")

    def decompress_file(self, src_path):
        dst_path = src_path[:-2] if src_path.endswith('.Z') else src_path[:-3]
        uc_func = unlzw.unlzw if src_path.endswith('.Z') else gzip.decompress
        with open(src_path, 'rb') as sf, open (dst_path, 'wb') as df:
            buffer = sf.read()
            df.write(uc_func(buffer))
        os.remove(src_path)
        return dst_path

    def crx2rnx(self, src_dir, src_path, doy, yy):
        sys = platform.system()
        if sys == "Windows":
            obs_name = src_path.split('\\')[-1]
            sta_name = str.lower(src_path.split('\\')[-1][0:4])
            re_name = sta_name + '{:03d}0.{:2d}o'.format(doy, yy)
            cmd = 'crx2rnx -f %s' % src_path
        else:
            obs_name = src_path.split('/')[-1]
            sta_name = str.lower(src_path.split('/')[-1][0:4])
            re_name = sta_name + '{:03d}0.{:2d}o'.format(doy, yy)
            cmd = 'CRX2RNX -f %s' % src_path
        print(cmd)
        subprocess.call(cmd, shell=True)
        rnx_path = src_dir + '/' + obs_name[:-3] + 'rnx'
        os.rename(rnx_path, os.path.join(src_dir, re_name))

        os.remove(src_path)

    def close(self):
        self.ftp.close()


def down_obs(ftp, down_sta_select, time):
    """
    down_sta_select: 若为[], 则默认下载多有mgex测站数据; 不为空且为四位的测站名则下载指定的测站数据；也可指定测站所属国，如CHN
    下载mgex的rinex3数据：ftp://igs.gnsswhu.cn/pub/gps/data/daily/
    长文件名改成短文件名 <SITE><DOY>.<YY>o
    数据存储在 /mgex/<YEAR>/<WEEK>/<DOY>/obs/
    """

    doy = common.ymd2doy(time.year, time.month, time.day)
    yy = time.year - 2000 if time.year >= 2000 else time.year - 1900
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)
    local_dir = LOCAL_DIR + str(time.year) + '/' + '{:d}/'.format(gps_week) \
                + '{:03d}/'.format(doy) + 'obs'
    remote_dir = 'pub/gps/data/daily/' + str(time.year) + '/' + '{:03d}'.format(doy) + '/' + '{:02d}d'.format(yy)
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    down_obs_files = ftp.file_classification(remote_dir, down_sta_select, 'obs')
    for file_name in down_obs_files:
        re_name = str.lower(file_name[:4]) + '{:03d}0.'.format(doy) + '{:2d}o'.format(yy)
        re_path = os.path.join(local_dir, re_name)
        if os.path.exists(re_path):
            print('%s is exits!' % re_path)
            continue
        local_path = os.path.join(local_dir, file_name)
        remote_path = os.path.join(remote_dir, file_name)
        ftp.down_file_single(local_path, remote_path)
        src_path = ftp.decompress_file(local_path)
        ftp.crx2rnx(local_dir, src_path, doy, yy)


def down_nav(ftp, time):
    """
    下载星历数据：ftp://igs.gnsswhu.cn/pub/gnss/mgex/daily/rinex3/<YEAR>/<DOY>/<YY>p/brdm<DOY>0.<YY>p.Z
    数据存储在 /mgex/<YEAR>/<WEEK>/<DOY>/prods/
    """
    doy = common.ymd2doy(time.year, time.month, time.day)
    yy = time.year - 2000 if time.year >= 2000 else time.year - 1900
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)

    local_dir = LOCAL_DIR + str(time.year) + '/' + '{:d}/'.format(gps_week) + '{:03d}/'.format(
        doy) + 'prods'
    remote_dir = 'pub/gnss/mgex/daily/rinex3/' + str(time.year) + '/' + '{:03d}'.format(
        doy) + '/' + '{:02d}p'.format(yy)
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    file_name = 'brdm%03d0.%02dp.Z' % (doy, yy)
    local_path = os.path.join(local_dir, file_name)
    if os.path.exists(local_path[:-2]):
        print("%s is exit!" % local_path[:-2])
        return
    remote_path = os.path.join(remote_dir, file_name)

    ftp.down_file_single(local_path, remote_path)
    ftp.decompress_file(local_path)


def down_igs_erpsnx(ftp, time):
    """
    下载igs周解的erp和snx文件: ftp://igs.gnsswhu.cn/pub/gps/products/<week>/igs<yy>P<week>.<snx, erp>.Z
    """
    doy = common.ymd2doy(time.year, time.month, time.day)
    yy = time.year - 2000 if time.year >= 2000 else time.year - 1900
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)
    local_dir = LOCAL_DIR + str(time.year) + '/' + '{:d}/'.format(gps_week)
    remote_dir = 'pub/gps/products/' + str(gps_week)
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    condition_list = list()
    condition_list.append('igs%dP%d.erp.Z' % (yy, gps_week))
    condition_list.append('igs%dP%d.snx.Z' % (yy, gps_week))

    down_erpsnx_files = ftp.file_classification(remote_dir, condition_list, 'igs')
    for file_name in down_erpsnx_files:
        local_path = os.path.join(local_dir, file_name)
        if os.path.exists(local_path[:-2]):
            print('%s is exits!' % local_path[:-2])
            continue
        remote_path = os.path.join(remote_dir, file_name)
        ftp.down_file_single(local_path, remote_path)
        ftp.decompress_file(local_path)


def down_pre(ftp, ac_list, time):
    """
    下载精密产品：ftp://igs.gnsswhu.cn/pub/gnss/products/mgex/<week>
    只下载MGEX的长文件名数据，即GPS 2038周 2019 1 28 后
    ac_list：发布机构 ['COD', 'GFZ', 'WUM']; 若为空则全部下载; GRG没提供erp暂不下载
    数据存储在 /mgex/<YEAR>/<DOY>/products/<ac>
    """
    doy = common.ymd2doy(time.year, time.month, time.day)
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)
    remote_dir = 'pub/gnss/products/mgex/' + str(gps_week)

    if ac_list[0] is 'GBM':
        remote_dir = 'GNSS/products/mgex/' + str(gps_week)
    if not ac_list:
        ac_list = ['COD', 'GFZ', 'WUM']

    for i in range(0, len(ac_list)):
        local_dir = LOCAL_DIR + str(time_idx.year) + '/' + '{:d}/'.format(gps_week) \
                    + str.lower(ac_list[i])
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)
        pre_files = ftp.file_classification(remote_dir, ac_list[i], 'pre')
        for file_name in pre_files:
            year = int(file_name[11:15])
            doy = int(file_name[15:18])
            gps_week, gps_dow = common.yrdoy2gpst(year, doy)
            re_path = rename_pre(file_name[:-3], ac_list[i], local_dir, gps_week, gps_dow)
            if os.path.exists(re_path):
                print('%s is exits!' % re_path)
                continue

            local_path = os.path.join(local_dir, file_name)
            remote_path = os.path.join(remote_dir, file_name)

            if ftp.down_file_single(local_path, remote_path) is False:
                continue
            dst_path = ftp.decompress_file(local_path)
            os.rename(dst_path, re_path)


def rename_pre(file_name, ac, src_dir, gps_week, gps_sow):
    if file_name.endswith('.SP3'):
        re_name = str.lower(ac) + str(gps_week) + str(gps_sow) + '.sp3'
    elif file_name.endswith('.CLK'):
        re_name = str.lower(ac) + str(gps_week) + str(gps_sow) + '.clk'
    elif file_name.endswith('.ERP'):
        re_name = str.lower(ac) + str(gps_week) + str(gps_sow) + '.erp'

    re_path = src_dir + '/' + re_name
    return re_path


def down_ion_dcb(ftp, time):
    doy = common.ymd2doy(time.year, time.month, time.day)
    yy = time.year - 2000 if time.year >= 2000 else time.year - 1900
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)

    condition_list = list()
    condition_list.append('CODG{:03d}0.{:02d}I.Z'.format(doy, yy))
    condition_list.append('P1C1{:2d}{:02d}.DCB.Z'.format(yy, time_idx.month))
    condition_list.append('P1P2{:2d}{:02d}.DCB.Z'.format(yy, time_idx.month))
    condition_list.append('P2C2{:2d}{:02d}_RINEX.DCB.Z'.format(yy, time_idx.month))
    local_dir = LOCAL_DIR + str(time.year) + '/' + str(gps_week) + '/' + '{:03d}'.format(doy) + '/' + 'prods'
    remote_dir = 'CODE/' + str(time.year)
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    ion_dcb_files = ftp.file_classification(remote_dir, condition_list, 'ion_dcb')
    for file_name in ion_dcb_files:
        if 'DCB' in file_name:
            local_dir = LOCAL_DIR + str(time.year) + '/' + 'dcb'
        else:
            local_dir = LOCAL_DIR + str(time.year) + '/' + str(gps_week) + '/' + '{:03d}'.format(
                doy) + '/' + 'prods'

        if not os.path.exists(local_dir):
            os.makedirs(local_dir)
        local_path = os.path.join(local_dir, file_name.split('/')[-1])
        if os.path.exists(local_path[:-2]):
            print("%s is exit!" % local_path[:-2])
            continue
        remote_path = remote_dir + '/' + file_name.split('/')[-1]
        ftp.down_file_single(local_path, remote_path)
        ftp.decompress_file(local_path)


def down_cas_dcb(ftp, time):
    doy = common.ymd2doy(time.year, time.month, time.day)
    yy = time.year - 2000 if time.year >= 2000 else time.year - 1900
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)

    local_dir = LOCAL_DIR + str(time.year) + '/' + '{:d}/'.format(gps_week) + '{:03d}/'.format(
        doy) + 'prods'
    remote_dir = 'product/dcb/mgex/' + str(time.year) + '/'
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    file_name = 'CAS0MGXRAP_%04d%03d0000_01D_01D_DCB.BSX.gz' % (time.year, doy)
    # file_name = 'CAS0MGXRAP_%04d%03d0000_01D_01D_OSB.BIA.gz' % (time.year, doy)
    local_path = os.path.join(local_dir, file_name)
    if os.path.exists(local_path[:-3]):
        print("%s is exit!" % local_path[:-3])
        return
    remote_path = os.path.join(remote_dir, file_name)

    ftp.down_file_single(local_path, remote_path)
    ftp.decompress_file(local_path)


def down_cod_osb(ftp, time):
    doy = common.ymd2doy(time.year, time.month, time.day)
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)
    remote_dir = 'pub/gnss/products/mgex/' + str(gps_week)
    local_dir = LOCAL_DIR + str(time.year) + '/' + '{:d}/'.format(gps_week) + '{:03d}/'.format(
        doy) + 'prods'
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    file_name = 'COD0MGXFIN_%04d%03d0000_01D_01D_OSB.BIA.gz' % (time.year, doy)
    local_path = os.path.join(local_dir, file_name)
    if os.path.exists(local_path[:-2]):
        print("%s is exit!" % local_path[:-2])
        return
    remote_path = os.path.join(remote_dir, file_name)
    ftp.down_file_single(local_path, remote_path)
    ftp.decompress_file(local_path)


def down_gbm_pre(ftp, time):
    doy = common.ymd2doy(time.year, time.month, time.day)
    gps_week, gps_dow = common.yrdoy2gpst(time_idx.year, doy)
    remote_dir = 'pub/gnss/products/mgex/' + str(gps_week)


DOWN_OBS_NAV_PRE = 1
DOWN_ION_DCB = 0
DOWN_CAS_DCB = 1
DOWN_GBM_PRE = 0
DOWN_COD_OSB = 0
sta_list = ['ABMF', 'DJIG', 'FAA1', 'HARB', 'JFNG', 'KARR', 'KIRU',
            'MAS1', 'MAYG', 'NIUM', 'NNOR', 'NYA2', 'PADO', 'SGOC',
            'SGPO', 'TOW2', 'URUM']

ac_list = ['COD', 'WUM', 'GRG']

if __name__ == "__main__":

    time_start = datetime.datetime(2020, 3, 15, 0, 0, 0)
    time_end = datetime.datetime(2020, 3, 15, 0, 0, 0)
    time_idx = time_start

    if DOWN_OBS_NAV_PRE:
        path = 'igs.gnsswhu.cn'
        ftp = FtpDown(path)
        ftp.login()

        while 1:
            if (time_idx - time_end).total_seconds() > 0:
                break

            down_obs(ftp, sta_list, time_idx)
            down_nav(ftp, time_idx)
            down_pre(ftp, ac_list, time_idx)
            down_igs_erpsnx(ftp, time_idx)
            time_idx = time_idx + relativedelta(days=1)
        ftp.close()

    if DOWN_ION_DCB:
        time_idx = time_start
        path = 'ftp.aiub.unibe.ch'
        ftp = FtpDown(path)
        ftp.login()

        while 1:
            if (time_idx - time_end).total_seconds() > 0:
                break

            down_ion_dcb(ftp, time_idx)
            time_idx = time_idx + relativedelta(days=1)
        ftp.close()

    if DOWN_CAS_DCB:
        time_idx = time_start
        path = 'ftp.gipp.org.cn'
        ftp = FtpDown(path)
        ftp.login()

        while 1:
            if (time_idx - time_end).total_seconds() > 0:
                break

            down_cas_dcb(ftp, time_idx)
            time_idx = time_idx + relativedelta(days=1)
        ftp.close()

    if DOWN_COD_OSB:
        time_idx = time_start
        path = 'igs.gnsswhu.cn'
        ftp = FtpDown(path)
        ftp.login()

        while 1:
            if (time_idx - time_end).total_seconds() > 0:
                break

            down_cod_osb(ftp, time_idx)
            time_idx = time_idx + relativedelta(days=1)
        ftp.close()

    if DOWN_GBM_PRE:
        time_idx = time_start
        path = 'ftp.gfz-potsdam.de'
        ftp = FtpDown(path)
        ftp.login()

        while 1:
            if (time_idx - time_end).total_seconds() > 0:
                break

            ac_list = ['GBM']
            down_pre(ftp, ac_list, time_idx)
            time_idx = time_idx + relativedelta(days=1)
        ftp.close()

    print("have fun")

