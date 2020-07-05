import math
import datetime
import os

SYS_GPS = 0x01
SYS_GLO = 0x04
SYS_GAL = 0x08
SYS_CMP = 0x20

MAX_GPS_SAT = 32
MAX_BDS_SAT = 45
MAX_GAL_SAT = 36
MAX_GLO_SAT = 27
MAX_QZS_SAT = 10
TOTAL_SAT = MAX_GPS_SAT + MAX_BDS_SAT + MAX_GAL_SAT + MAX_GLO_SAT + MAX_BDS_SAT

FONT_SIZE = 10
# 黑红搭配, 经典万能，红色用于重点数据线，黑色用于参考线
COLOR_RED = (255, 59, 59)
COLOR_BLACK = (7, 7, 7)

# 红蓝配色， 看起来不会很死板
COLOR_RED1 = (254, 129, 125)
COLOR_BLUE = (129, 184, 223)

# 经典三色搭配
COLOR3_1 = (77, 133, 189)
COLOR3_2 = (247, 144, 61)
COLOR3_3 = (89, 169, 90)

# 经典三色搭配，一暖二冷突出重点
COLOR3_11 = (210, 32, 39)
COLOR3_21 = (56, 89, 137)
COLOR3_31 = (127, 165, 183)


# 四色搭配
COLOR4_1 = (23, 23, 23)
COLOR4_2 = (6, 223, 6)
COLOR4_3 = (255, 28, 0)
COLOR4_4 = (0, 37, 255)

# 五色搭配应尽量降低饱和度
COLOR5_1 = (1, 86, 153)
COLOR5_2 = (250, 192, 15)
COLOR5_3 = (243, 118, 74)
COLOR5_4 = (95, 198, 201)
COLOR5_5 = (79, 89, 109)

# 六色搭配
COLOR6_1 = (203/256, 180/256, 123/256)
COLOR6_2 = (91/256, 183/256, 205/256)
COLOR6_3 = (71/256, 120/256, 185/256)
COLOR6_4 = (84/256, 172/256, 117/256)
COLOR6_5 = (197/256, 86/256, 89/256)
COLOR6_6 = (117/256, 114/256, 181/256)


class gtime:
    def __init__(self):
        gtime.time = 0.0
        gtime.sec = 0.0


def get_filelist(path):
    file_list = []
    for home, dirs, files in os.walk(path):
        for file_name in files:
            file_list.append(os.path.join(home, file_name))
    return file_list


# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
def date_diff_seconds(date):
    date_1970jan1 = datetime.datetime(1970, 1, 1, 0, 0, 0)
    time_delta = date - date_1970jan1
    return time_delta.days*24*3600 + time_delta.seconds


# convert 2-digit year to 4-digit year
def yy2yyyy(yy):
    if yy >= 1900:
        yyyy = yy
    elif 50 < yy < 1900:
        yyyy = yy + 1900
    elif yy <= 50:
        yyyy = yy + 2000
    return yyyy


# convert year, day of year to GPS week, day of week
def yrdoy2gpst(year, doy):
    year0 = yy2yyyy(year)
    date_1980jan6 = datetime.datetime(1980, 1, 6, 0, 0, 0)
    date = datetime.datetime(year0, 1, 1, 0, 0, 0)
    time_delta = date - date_1980jan6
    days_delta = time_delta.days + doy - 1
    gps_week = int(days_delta/7)
    gps_dow = int(days_delta - gps_week*7)
    return gps_week, gps_dow


# convert GPS week, seconds of week to modified Julian date(MJD)
def gpst2mjd(gps_week, gps_sow):
    mjd = gps_week*7 + 44244 + gps_sow/86400
    return mjd


# convert year, day of year to GPS week, day of week
def gpst2yrdoy(gps_week, gps_sow):
    mjd = gpst2mjd(gps_week, gps_sow)
    year, doy = mjd2yrdoy(mjd)
    return year, doy


def ymd2doy(year, mon, day):
    dn = datetime.datetime(year, mon, day, 0, 0, 0)
    return int(dn.strftime("%j"))


# convert year, month, day or year, day of year to modified Julian date(MJD)
def ymd2mjd(year, month, day):
    doy_of_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

    year0 = yy2yyyy(year)
    if year0 < 0 or month < 0 or month > 12 or day > 336 or (month != 0 and day > 31):
        print('*** ERROR(ymd2mjd): Incorrect date(year, month, day):', year0, month, day)

    if month == 0:
        im, id0 = yrdoy2ymd(year0, day)
    else:
        im = month
        id0 = day
    year1 = year0
    if im <= 2:
        year1 = year1 - 1
    mjd = 365 * year0 - 678941 + int(year1 / 4) - int(year1 / 100) + int(year1 / 400) + id0
    im = im - 1
    if im != -1:
        mjd = mjd + doy_of_month[im]
    return mjd


# convert year, day of year to month, day
def yrdoy2ymd(year, doy):
    year0 = yy2yyyy(year)
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (year0 % 4 == 0 and year0 % 100 != 0) or year0 % 400 == 0:
        days_in_month[1] = 29
    id0 = doy
    for imonth in range(0, 12):
        id0 = id0 - days_in_month[imonth]
        if id0 > 0:
            continue
        iday = id0 + days_in_month[imonth]
        break
    imonth += 1
    return imonth, iday


#  convert MJD to year, day of year
def mjd2yrdoy(mjd):
    year = int((mjd + 678940) / 365)
    doy = mjd - ymd2mjd(year, 1, 1)
    while doy <= 0:
        year = year - 1
        doy = mjd - ymd2mjd(year, 1, 1) + 1
    return year, doy


# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
def get_sat_id(sat_prn):

    if 0 <= sat_prn <= 31:
        if 0 <= sat_prn <= 8:
            sat_id = 'G' + '0' + str(sat_prn + 1)
        else:
            sat_id = 'G' + str(sat_prn + 1)
    elif 32 <= sat_prn <= 76:
        if 0 <= sat_prn - 32 <= 8:
            sat_id = 'C' + '0' + str(sat_prn - MAX_GPS_SAT + 1)
        else:
            sat_id = 'C' + str(sat_prn - MAX_GPS_SAT + 1)
    elif 77 <= sat_prn <= 112:
        if 0 <= sat_prn - 77 <= 8:
            sat_id = 'E' + '0' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT + 1)
        else:
            sat_id = 'E' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT + 1)
    elif 113 <= sat_prn <= 139:
        if 0 <= sat_prn - 113 <= 8:
            sat_id = 'R' + '0' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT - MAX_GAL_SAT + 1)
        else:
            sat_id = 'R' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT - MAX_GAL_SAT + 1)
    elif 140 <= sat_prn <= 149:
        if 0 <= sat_prn - 140 <= 8:
            sat_id = 'J' + '0' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT - MAX_GAL_SAT - MAX_GLO_SAT + 1)
        else:
            sat_id = 'J' + str(sat_prn - MAX_GPS_SAT - MAX_BDS_SAT - MAX_GAL_SAT - MAX_GLO_SAT + 1)
    return sat_id


def get_max_set_sat(sys):
    sat_num = MAX_GPS_SAT
    if sys & SYS_GLO:
        sat_num += MAX_GLO_SAT
    if sys & SYS_GAL:
        sat_num += MAX_GAL_SAT
    if sys & SYS_CMP:
        sat_num += MAX_BDS_SAT

    return sat_num


# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
def get_sat_no(sat_id):
    sat_sys = sat_id[0]
    sat_no = int(sat_id[1:3])
    if sat_sys == 'C':
        sat_no += MAX_GPS_SAT
    elif sat_sys == 'E':
        sat_no += MAX_GPS_SAT + MAX_BDS_SAT
    elif sat_sys == 'R':
        sat_no += MAX_GPS_SAT + MAX_BDS_SAT + MAX_GAL_SAT
    elif sat_sys == 'J':
        sat_no += MAX_GPS_SAT + MAX_BDS_SAT + MAX_GAL_SAT + MAX_GLO_SAT

    return sat_no


def get_average(records):
    """
    平均值
    """
    return sum(records)/len(records)


def get_variance(records):
    """
    方差: 反映一个数据的离散程度
    """
    average = get_average(records)
    return sum([(x-average) ** 2 for x in records])/len(records)


def get_standard_deviation(records):
    """
    标准差: ==均方差 反映一个数据集的离散程度
    """
    variance = get_variance(records)
    return math.sqrt(variance)


def get_rms(records):
    """
    均方根值：反映有效值而不是平均值, 不一定反映的是误差，可以是一个序列
    均方误差
    """
    return math.sqrt(sum([x ** 2 for x in records])/len(records))


def get_mse(records):
    """
    均方误差：估计值与真值偏差
    """
    return sum([x ** 2 for x in records])/len(records)


def get_rmse(records):
    """
    均方根误差： 是均方根误差的算术平方根
    """
    mse = get_mse(records)
    if mse:
        return math.sqrt(mse)
    else:
        return None