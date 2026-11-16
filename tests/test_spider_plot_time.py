import sys,unittest
from pathlib import Path
from datetime import datetime,timezone
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spider_plot_time import date_coordinates,format_time_axis

class TimeAxisTests(unittest.TestCase):
    def test_local_midnight_rollover(self):
        stamps=[datetime(2018,8,21,21,30,tzinfo=timezone.utc),datetime(2018,8,21,22,30,tzinfo=timezone.utc)]
        x=date_coordinates(stamps);fig,ax=plt.subplots();format_time_axis(ax,'Europe/Stockholm')
        fmt=ax.xaxis.get_major_formatter()
        self.assertEqual(fmt(x[0]),'21 Aug\n2018\n23:30')
        self.assertEqual(fmt(x[1]),'22 Aug\n2018\n00:30')
        self.assertAlmostEqual(x[1]-x[0],1/24)
        plt.close(fig)
    def test_reject_naive_timestamps(self):
        with self.assertRaises(ValueError):date_coordinates([datetime(2018,8,21)])

if __name__=='__main__':unittest.main()
