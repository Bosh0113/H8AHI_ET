from datetime import datetime, timedelta
import os
# import subprocess

START_TIME = '2018-07-19T00:00:00Z' # local time
END_TIME = '2018-07-19T23:59:59Z'

UTC_OFFSET = 9 # hour
time_internal = 10  # mins

storage_path = r'D:\PhD_Workspace\H8_ET\231017\JP_AMATERASS'


def download_AMATERASS_data(file_time):

    file_suffixes = ['.dwn.sw.flx.sfc.msm.1km.bin.bz2', '.rh.sfc.msm.1km.bin.bz2', '.tsfc.msm.1km.bin.bz2'] # SRd, RH, Ta
    
    clear_flag = 0
    for file_suffix in file_suffixes:
        file_time_str = file_time
        file_name = file_time_str + file_suffix
        local_file = storage_path + '/' + file_name

        if not os.path.exists(local_file):
            clear_flag = 1
    if clear_flag:
        for file_suffix in file_suffixes:
            file_time_str = file_time
            file_name = file_time_str + file_suffix
            local_file = storage_path + '/' + file_name

            if os.path.exists(local_file):
                os.remove(local_file)
                print(local_file)


if __name__ == "__main__":

    data_start_date = datetime.strptime(START_TIME, "%Y-%m-%dT%H:%M:%SZ") - timedelta(hours=UTC_OFFSET)
    data_end_date = datetime.strptime(END_TIME, "%Y-%m-%dT%H:%M:%SZ") - timedelta(hours=UTC_OFFSET)

    temp_data = data_start_date
    while temp_data < data_end_date:
        current_time_str = temp_data.strftime("%Y%m%d%H%M")
        download_AMATERASS_data(current_time_str)
        temp_data = temp_data + timedelta(minutes=time_internal)