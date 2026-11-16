"""Shared dated local-time axes for spider night plots; computations stay in UTC."""
from zoneinfo import ZoneInfo
import matplotlib.dates as mdates


def date_coordinates(stamps):
    if any(t.tzinfo is None or t.utcoffset() is None for t in stamps):raise ValueError('Timezone-aware timestamps required')
    return mdates.date2num(stamps)


def format_time_axis(ax,zone):
    zone=ZoneInfo(zone) if isinstance(zone,str) else zone
    locator=mdates.AutoDateLocator(tz=zone,minticks=3,maxticks=4)
    locator.intervald[mdates.HOURLY]=[1,2,3,4,6,12]
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b\n%Y\n%H:%M',tz=zone))
    ax.set_xlabel(f'Local time ({zone.key})',fontsize=9)
    ax.tick_params(axis='x',labelsize=8)
