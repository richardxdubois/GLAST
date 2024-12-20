import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from bokeh.models.widgets import Div, NumberFormatter
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import ColumnDataSource, Span, LinearAxis, Range1d, LinearColorMapper, Whisker
from bokeh.plotting import figure, output_file, reset_output, show, save
from bokeh.layouts import row, layout, column

# Define the time range
time_days = 3650  # Example: 10 years
time = np.arange(0, time_days, 1)  # Daily steps

# Define parameters for the modulations
base_value = 10  # Constant base value
amplitude_1 = 3  # Amplitude for the first period (26.496 days)
period_1 = 26.496  # First modulation period in days

amplitude_2 = 1  # Amplitude for the second period (1667 days)
period_2 = 1667  # Second modulation period in days

# Generate the time series
modulation_1 = amplitude_1 * np.sin(2 * np.pi * time / period_1)
modulation_2 = amplitude_2 * np.sin(2 * np.pi * time / period_2)
time_series = base_value + modulation_1 + modulation_2

# Create a DataFrame for better handling and visualization
time_series_df = pd.DataFrame({
    'Time (days)': time,
    'Value': time_series
})

# Plot the time series
p = figure(title="Time series", x_axis_label='Time (days)', y_axis_label='amplitude', width=750)
p.line(x=time_series_df['Time (days)'], y=time_series_df['Value'])
p.y_range.start = 0

frequency = np.linspace(1/28., 1/25., 1000)  # for orbital 1000
freq_days = 1. / frequency

power = LombScargle(t=time_series_df['Time (days)'], y=time_series_df['Value']).power(frequency)
peaks_ls, props_ls = find_peaks(power, height=0.1 * max(power))
pk_days = (1. / frequency[peaks_ls])

f1 = figure(title="full time span: power vs frequency",
            x_axis_label='period (days)', y_axis_label='power',
            width=750)
f1.line(freq_days, power, line_width=2)

lay = layout(p, f1)
save(lay, title="simulation orb & super")
