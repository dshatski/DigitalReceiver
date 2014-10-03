#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: USRP1 digital receiver
# Author: Dmitry Shatskiy
# Description: USRP1 to file
# Generated: Thu Apr 26 23:37:13 2012
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from optparse import OptionParser
import time
import os

class usrp1(gr.top_block):

	def __init__(self, f0_lo=149.97e6, f0_hi=399.92e6, path150="/home/alfagot/Desktop/Passes/_data150.dat", path400="/home/alfagot/Desktop/Passes/_data400.dat", t0=660):
		gr.top_block.__init__(self, "USRP1 digital receiver")

		##################################################
		# Parameters
		##################################################
		self.f0_lo = f0_lo
		self.f0_hi = f0_hi
		self.path150 = path150
		self.path400 = path400
		self.t0 = t0

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 250e3
		self.decim = decim = 2

		##################################################
		# Blocks
		##################################################
		self.uhd_usrp_source_0 = uhd.usrp_source(
			device_addr="type=usrp1",
			stream_args=uhd.stream_args(
				cpu_format="fc32",
				channels=range(2),
			),
		)
		self.uhd_usrp_source_0.set_subdev_spec("B:RX1 B:RX2", 0)
		self.uhd_usrp_source_0.set_samp_rate(samp_rate)
		self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(f0_lo, 0), 0)
		self.uhd_usrp_source_0.set_gain(0, 0)
		self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(f0_hi, 0), 1)
		self.uhd_usrp_source_0.set_gain(0, 1)
		self.low_pass_filter_1 = gr.fir_filter_ccf(decim, firdes.low_pass(
			1, samp_rate, 12e3, 1e3, firdes.WIN_HAMMING, 6.76))
		self.low_pass_filter_0 = gr.fir_filter_ccf(decim, firdes.low_pass(
			1, samp_rate, 5e3, 1e3, firdes.WIN_HAMMING, 6.76))
		self.gr_file_sink_1 = gr.file_sink(gr.sizeof_gr_complex*1, path400)
		self.gr_file_sink_1.set_unbuffered(False)
		self.gr_file_sink_0 = gr.file_sink(gr.sizeof_gr_complex*1, path150)
		self.gr_file_sink_0.set_unbuffered(False)

		##################################################
		# Connections
		##################################################
		self.connect((self.uhd_usrp_source_0, 0), (self.low_pass_filter_0, 0))
		self.connect((self.low_pass_filter_0, 0), (self.gr_file_sink_0, 0))
		self.connect((self.low_pass_filter_1, 0), (self.gr_file_sink_1, 0))
		self.connect((self.uhd_usrp_source_0, 1), (self.low_pass_filter_1, 0))

	def get_f0_lo(self):
		return self.f0_lo

	def set_f0_lo(self, f0_lo):
		self.f0_lo = f0_lo
		self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(self.f0_lo, 0), 0)

	def get_f0_hi(self):
		return self.f0_hi

	def set_f0_hi(self, f0_hi):
		self.f0_hi = f0_hi
		self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(self.f0_hi, 0), 1)

	def get_path150(self):
		return self.path150

	def set_path150(self, path150):
		self.path150 = path150

	def get_path400(self):
		return self.path400

	def set_path400(self, path400):
		self.path400 = path400

	def get_t0(self):
		return self.t0

	def set_t0(self, t0):
		self.t0 = t0

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.low_pass_filter_1.set_taps(firdes.low_pass(1, self.samp_rate, 12e3, 1e3, firdes.WIN_HAMMING, 6.76))
		self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
		self.low_pass_filter_0.set_taps(firdes.low_pass(1, self.samp_rate, 5e3, 1e3, firdes.WIN_HAMMING, 6.76))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--f0-lo", dest="f0_lo", type="eng_float", default=eng_notation.num_to_str(149.97e6),
		help="Set f0_lo [default=%default]")
	parser.add_option("", "--f0-hi", dest="f0_hi", type="eng_float", default=eng_notation.num_to_str(399.92e6),
		help="Set f0_hi [default=%default]")
	parser.add_option("", "--path", dest="path", type="string", default="/home/alfagot/Desktop/Passes",
		help="Set path [default=%default]")
	parser.add_option("", "--t0", dest="t0", type="intx", default=660,
		help="Set t0 [default=%default]")
	(options, args) = parser.parse_args()
	try:
	    os.makedirs(options.path)
	except OSError:
	    pass
	tb = usrp1(f0_lo=options.f0_lo, f0_hi=options.f0_hi, path150=options.path+"/_data150.dat", path400=options.path+"/_data400.dat", t0=options.t0)
	tb.start()
	start_time = time.time()  # remember when we've started
	while (time.time() - start_time) < options.t0:
		pass # do nothing
	tb.stop()
