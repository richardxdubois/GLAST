from math import pi
import numpy as np

import pandas as pd
import pickle

from bokeh.models import BasicTicker, PrintfTickFormatter, ColumnDataSource
from bokeh.plotting import figure, show
from bokeh.transform import linear_cmap
#from bokeh.palettes import Viridis256 as palette
from bokeh.palettes import Plasma256 as palette

with open('/Users/richarddubois/Code/GLAST/tmp/heatmap.pkl', 'rb') as f:
    data = pickle.load(f)

all_x = data[0]
all_y = data[1]
all_flux = data[2]

p_bins = np.arange(10)
q_bins = p_bins
q_bins = np.append(q_bins, p_bins)
dict_ticker = {}
for i in all_x:
    dict_ticker[all_x[i]] = str(q_bins[i]/10.)

source = ColumnDataSource(data=dict(x=all_x, y=all_y, flux=all_flux))

# this is the colormap from the original NYTimes plot
colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]

TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

p = figure(title="Orbital phase calculated per super phase bin - flux vs super and orbital phases",
           x_axis_location="above", width=900, height=600,
           tools=TOOLS, toolbar_location='below', y_axis_label="super phase", x_axis_label="orbital phase",
           tooltips=[('phases', 'Orbital: @x Super: @y'), ('flux', '@flux')])

p.grid.grid_line_color = None
p.axis.axis_line_color = None
p.axis.major_tick_line_color = None
p.axis.major_label_text_font_size = "12px"
p.axis.major_label_standoff = 0
p.xaxis.major_label_orientation = pi / 3

r = p.rect(x="x", y="y", width=1, height=1, source=source,
           fill_color=linear_cmap("flux", palette=palette, low=0., high=1.2e-6),
           line_color=None, )
p.xaxis.ticker = all_x
p.xaxis.major_label_overrides = dict_ticker
p.yaxis.ticker = all_y
p.yaxis.major_label_overrides = dict_ticker

p.add_layout(r.construct_color_bar(
    major_label_text_font_size="12px",
    ticker=BasicTicker(desired_num_ticks=len(colors)),
    #formatter=PrintfTickFormatter(format="%d%%"),
    label_standoff=6,
    border_line_color=None,
    padding=5,
), 'right')

show(p)
