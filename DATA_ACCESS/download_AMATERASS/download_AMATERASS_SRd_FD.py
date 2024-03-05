from datetime import datetime, timedelta
from ftplib import FTP
import os
# import subprocess

START_TIME = '2018-01-01T00:00:00Z' # local time
END_TIME = '2019-12-31T23:59:59Z'

UTC_OFFSET = 9 # hour
time_internal = 10  # mins

storage_path = os.getcwd()


def download_AMATERASS_data(file_time):

    file_suffixes = ['.dwn.sw.flx.sfc.fld.4km.bin.bz2'] # SRd
    
    for file_suffix in file_suffixes:
        file_time_str = file_time
        file_name = file_time_str + file_suffix
        ftp_path= '/quasi-realtime/himawari829/archived/FD/'+ file_time_str[:6] + '/' + file_time_str[:8] + '/' + file_name
        remote_url = 'ftp://amaterass.cr.chiba-u.ac.jp' + ftp_path

        storage_folder = os.path.join(storage_path, file_time_str[:6], file_time_str[:8])
        if not os.path.exists(storage_folder):
            os.makedirs(storage_folder)

        local_file = storage_folder + '/' + file_name

        if not os.path.exists(local_file):
            try:
                with open(local_file, 'wb') as f:
                    ftp.retrbinary('RETR ' + ftp_path, f.write, 1024 * 1024)
                # p = subprocess.Popen('lbzip2 -d {}'.format(local_file), shell=True)
                # p.communicate()
            except Exception as e:
                print(remote_url)
                print(e)
                os.remove(local_file)


if __name__ == "__main__":

    ftp = FTP()
    ftp.connect('amaterass.cr.chiba-u.ac.jp', 21)
    ftp.login()

    data_start_date = datetime.strptime(START_TIME, "%Y-%m-%dT%H:%M:%SZ") - timedelta(hours=UTC_OFFSET)
    data_end_date = datetime.strptime(END_TIME, "%Y-%m-%dT%H:%M:%SZ") - timedelta(hours=UTC_OFFSET)

    temp_data = data_start_date
    while temp_data < data_end_date:
        current_time_str = temp_data.strftime("%Y%m%d%H%M")
        download_AMATERASS_data(current_time_str)
        temp_data = temp_data + timedelta(minutes=time_internal)

    ftp.close()