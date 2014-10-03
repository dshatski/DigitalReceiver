#!/usr/bin/python
import sys, os, datetime
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import *
from scipy import io, optimize

def peak(i):
    y = 8485.55446-1073.1255*tanh(0.000449676974*(i-103.744525))
    return int(y)

def TECanalysis():
    A_coef = 80.6
    c_speed = 2.998e8   # speed of light
    fr = 50.0e6         # basic frequency
    R=6378.14 #Earth radius, km
    hi=300 #height of ionosphere in a singe thin layer model, km
    h2407=989.2 #average altitude of Cosmos 2407, km
    h2414=948.5 #average altitude of Cosmos 2414, km
    h2429=993.3 #average altitude of Cosmos 2429, km
    T2407=104.7*60 #orbital period of Cosmos 2407, s
    T2414=103.8*60 #orbital period of Cosmos 2414, s
    T2429=104.7*60 #orbital period of Cosmos 2429, s
    ratio2407=14997.0/39992.0 #central frequencies ratio for Cosmos 2407
    ratio2414=14997.0/39992.0 #central frequencies ratio for Cosmos 2414
    ratio2429=15003.0/40008.0 #central frequencies ratio for Cosmos 2429
    halfband=16#halfband for the 400MHz ifft (should be divisible by 8)
    file150 = '_data150.dat'
    file400 = '_data400.dat'
    if os.path.exists(file400)& os.path.exists(file150):
        data_size = os.stat(file400)[6] # Get filesisze
        data_size = min(data_size, os.stat(file150)[6]) # Get filesize
    else:
        return
    sample_f = 125000.0
    nfft = 16384
    delta_f = sample_f / float(nfft)  #frequency resolution
#########################################################################################
    skip_sec = 0.0
    read_sec = 999.0
    maxelev = 80.0	#max elevation angle, deg
    sat = 2414
    if sat==2407:
	hs=h2407
	Ts=T2407
	ratio=ratio2407
    if sat==2414:
	hs=h2414
	Ts=T2414
	ratio=ratio2414
    if sat==2429:
	hs=h2429
	Ts=T2429
	ratio=ratio2429
    ws=2*pi/Ts #angular velocity, rad/s
    time0 = 734.0/2-skip_sec
    theta0 = -time0*ws #starting zenith angle in equatorial coordinates, rad
    print 'Starting zenith angle:',theta0/pi*180, 'deg. Initializing arrays...'
#########################################################################################
    if read_sec == 999.0:
        read_sec = float((data_size // (8 * sample_f)) - skip_sec)
    nskip_bytes = int((sample_f * skip_sec) // nfft) * nfft * 8
    ntdata = int((sample_f * read_sec) // nfft) * nfft
    nparam = ntdata / nfft * 2 -1 # Number of spectral data (half interleave)
    npfact = 128
    npst = nfft/npfact/4
    n_phase = ntdata / npfact      #Number of phase data for TEC evaluation (125000*128*<length_in_seconds>//16384)
    freq = arange(nfft, dtype = 'f') * delta_f - sample_f/2.0
    buf = arange(ntdata, dtype = 'F')
    ser = arange(nfft, dtype = 'D')
    fft_result = arange(nfft, dtype = 'D')
    power = arange(nfft, dtype = 'd')
    if150 = arange(nparam, dtype = 'l')
    if400 = arange(nparam, dtype = 'l')
    if400_a = arange(nparam, dtype = 'l')
    if400_b = arange(nparam, dtype = 'l')
    freq400 = arange(nparam, dtype = 'f')
    sigpow150 = arange(nparam, dtype = 'f')
    sigpow400 = arange(nparam, dtype = 'f')
    fil_ser150 = arange(n_phase, dtype = 'F')
    fil_ser400 = arange(n_phase, dtype = 'F')
    dif_p = arange(n_phase, dtype = 'f')
    pTEC = arange(n_phase, dtype = 'f')
    differ = arange(n_phase, dtype = 'f') * 0.0 # store the derivative
    phase_sec = arange(n_phase, dtype = 'f') * (float(npfact)/sample_f) + skip_sec
    param_sec = linspace(skip_sec, (nparam-1)*(float(0.5*nfft*ntdata/(ntdata-0.5*nfft))/sample_f) + skip_sec, nparam)
    theta = abs((phase_sec-skip_sec) * ws + theta0)
    d=sqrt(R*R+(R+hs)*(R+hs)-2*R*(R+hs)*cos(theta))
    slant_phi=arccos((R+hs)*sin(theta)/d)
    phi=arcsin(sin(slant_phi)*sin(maxelev/180.0*pi)) #take into account the satellite orbit inclination
    ksi=arcsin(R/(R+hi)*cos(phi))
    cos_ksi=cos(ksi)
    sec=1/cos_ksi
    arg=sec.argmin()

    plot(phase_sec,d)
    xlabel('time (s)')
    ylabel('distance (km)')
    title('distance to the satellite')
    savefig('_dist.png')
    close()

    plot(phase_sec,phi/pi*180)
    xlabel('time (s)')
    ylabel('phi (deg)')
    title('elevation angle')
    savefig('_phi.png')
    close()

    plot(phase_sec,sec)
    xlabel('time (s)')
    ylabel('sec(ksi)')
    title('secant effect')
    savefig('_sec.png')
    close()

    #print 'elevation at the start:',(slant_phi[0]/pi*180),' deg (',(phi[0]/pi*180),' deg after correcting for inclination)'
    #print 'elevation at the end:',(slant_phi[n_phase-1]/pi*180),' deg (',(phi[n_phase-1]/pi*180),' deg after correcting for inclination)'

    print 'Skipping',int(skip_sec),'seconds (',nskip_bytes,'bytes), reading',int(read_sec),'seconds containing',int(sample_f * read_sec),'samples, i.e.',int((sample_f * read_sec) // nfft),'bins of',nfft,'samples each for a total of',ntdata,'samples'

#########################################################################################
    f400 = open(file400, 'rb')
    f400.seek(nskip_bytes)
    print 'Buffering 400MHz data..'
    buf = io.fread(f400, ntdata, 'F')
    for ip in range(0, nparam):
        noffset = ip * nfft / 2
	if (ip % 2000) == 0:
		print int(float(ip)/nparam*25),'% complete'
        ser = buf[noffset: noffset+nfft]
	ser=ser*nfft*0.5 # scale the data so that the max value is 1
        fft_result = fft(ser)/nfft
        fft_result = fftpack.fftshift(fft_result)   #Rotate freq.
	fft_result[8192] = 0.0 + 0.0j
        power = fft_result.real **2 + fft_result.imag **2     #Signal power
	if (ip == 0):
		if400[ip] = power.argmax()
	if (ip > 0):
		if_init=if400[ip-1]-halfband #start one halfband off of the previous result
		if400[ip] = power[if_init: if_init+(2*halfband+1)].argmax() + if_init #look around +-halfband of the previous result#8514
    del buf
    f400.close
#########################################################################################
    fitfunc = lambda p, x: p[3] - p[0] * tanh(p[1] * (x - p[2])) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

    P0 = [1070., 0.007, phase_sec[n_phase/2], 8216] # Initial guess for the parameters

    P1, success = optimize.leastsq(errfunc, P0[:], args=(param_sec, if400))
    #P1 = [9.11939377e+02, 4.50602823e-03, 3.14300432e+02, 8.49913297e+03]
    if400_a = fitfunc(P1, param_sec)
    #print P1

    for ip in range(0, nparam):
	if (if400[ip-1]-if400[ip])>1:
	    if400_b[ip]=int(if400_a[ip])
	else:
	    if400_b[ip]=if400[ip]

    figure()
    plot(param_sec, abs(if400_b-if400), 'r')
    savefig('_dif.png')
    close()

    figure()
    plot(param_sec, if400) # Plot the data and the fit
    plot(param_sec, if400_b) # Plot the data and the fit
    savefig('_initguess.png')
    close()

    for ip in range(0, nparam):
	if400[ip]=int(if400_b[ip])
#########################################################################################
    f400 = open(file400, 'rb')
    f400.seek(nskip_bytes)
    buf = io.fread(f400, ntdata, 'F')
    fil_ser400[0:npst] = 1.0 + 0.0j
    fil_ser400[-npst:] = 1.0 + 0.0j
    ifst = nfft/4
    ifen = ifst * 3
    for ip in range(0, nparam):
        noffset = ip * nfft / 2
	if (ip % 2000) == 0:
		print 25+int(float(ip)/nparam*25),'% complete'
        ser = buf[noffset: noffset+nfft]
	ser=ser*nfft*0.5 # scale the data so that the max value is 1
        fft_result = fft(ser)/nfft
        fft_result = fftpack.fftshift(fft_result)   #Rotate freq.
        power = fft_result.real **2 + fft_result.imag **2     #Signal power
        sigpow400[ip] = 0.0
        freq400[ip] = 0.0
        for j in range(if400[ip]-halfband, if400[ip]+(halfband+1)):#        for j in range(if400[ip]-2, if400[ip]+3):
            sigpow400[ip] = sigpow400[ip] + power[j]
            freq400[ip] = freq400[ip] + freq[j] * power[j]
        freq400[ip] = freq400[ip] / sigpow400[ip]
        ist = ip*npst*2 + npst
        ien = ist + npst*2
        if if400[ip]<8176 or if400[ip]>8208
            fft_result[:if400[ip]-halfband] = 0.0 + 0.0j
            fft_result[if400[ip]+(halfband+1):] = 0.0 + 0.0j
            fft_result = fftpack.ifftshift(fft_result)  #Rotate freq. back
            ser = ifft(fft_result)
            fil_ser400[ist:ien] = ser[ifst:ifen:npfact]
        else:
	    print 'Filling with dummy signal'
            fil_ser400[ist:ien] = 1.0 + 0.0j
    del buf
#########################################################################################
    f150 = open(file150, 'rb')
    f150.seek(nskip_bytes)
    print 'Buffering 150MHz data..'
    f150.seek(nskip_bytes, 0)
    buf = io.fread(f150, ntdata, 'F')
    fil_ser150[0:npst] = 1.0 + 0.0j
    fil_ser150[-npst:] = 1.0 + 0.0j
    ifst = nfft/4
    ifen = ifst * 3
    halfband=int(halfband/8*3) #set the halfband for 150MHz
    for ip in range(0, nparam):
        noffset = ip * nfft / 2
	if (ip % 1000) == 0:
		print 50+int(float(ip)/nparam*50),'% complete'

        ser = buf[noffset: noffset+nfft]
	ser=ser*nfft*0.5 # scale the data so that the max value is 1
        fft_result = fft(ser)/nfft
        fft_result = fftpack.fftshift(fft_result)   #Rotate freq.
        power = fft_result.real **2 + fft_result.imag **2     #Signal power
        if isnan(abs(freq400[ip])):
	     print 'freq400[',ip,'] = NaN, excluding...'
             freq400[ip]=freq400[ip-1] #exclude the NaNs
        if_init = int(((freq400[ip]*ratio) + sample_f/2.0)/delta_f) - halfband
        if150[ip] = power[if_init: if_init+2*halfband+1].argmax() + if_init
        sigpow150[ip] = 0.0
        for j in range(if150[ip]-halfband, if150[ip]+(halfband+1)):
            sigpow150[ip] = sigpow150[ip] + power[j]
        ist = ip*npst*2 + npst
        ien = ist + npst*2
        if if150[ip]<8186 or if150[ip]>8198: 
            fft_result[:if150[ip]-halfband] = 0.0 + 0.0j
            fft_result[if150[ip]+(halfband+1):] = 0.0 + 0.0j
            fft_result = fftpack.ifftshift(fft_result)  #Rotate freq. back
            ser = ifft(fft_result)
            fil_ser150[ist:ien] = ser[ifst:ifen:npfact]
        else:
	    print 'Filling with dummy signal'
            fil_ser150[ist:ien] = 1.0 + 0.0j
    del buf
#########################################################################################
    print 'Processing data..'
    for ip in range(0, n_phase):
        dif_p[ip] = angle((fil_ser150[ip] ** 8) / (fil_ser400[ip] ** 3))
    pTEC = unwrap(dif_p)#pTEC = unwrap(dif_p)/8.0
    coef = pi * A_coef / fr / c_speed * 55.0 / 24.0
    pTEC = pTEC / coef
    f_out1 = open('_TEC1.txt', 'w')
    f_out1.write('Non-corrected non-calibrated slant TEC (Total %d lines)\n' % n_phase)
    for ip in range(0, n_phase):
        write_buf = '%g %g %g %g\n' % (phase_sec[ip], pTEC[ip], fil_ser150[ip], fil_ser400[ip])
        f_out1.write(write_buf)
    f_out1.close

    span=min(arg, 10000) #if overhead point is closer than 10000 samples
    (correction,b) = polyfit(phase_sec[arg-span:arg+span],pTEC[arg-span:arg+span],1)
    print 'correction =', correction/1.0e16, 'TECU/s'
    #plot the overhead region
    plot(phase_sec[arg-span:arg+span],pTEC[arg-span:arg+span]/1.0e16)
    plot(phase_sec[arg-span:arg+span],(phase_sec[arg-span:arg+span]*correction+b)/1.0e16)
    xlabel('time (s)')
    ylabel('relative TEC (TECU)')
    title('Overhead slope fitting')
    savefig('_corr.png')
    close()
 
    #plot the measured TEC
    plot(phase_sec,pTEC/1.0e16)
    plot(phase_sec,(phase_sec*correction+b)/1.0e16)
    xlabel('time (s)')
    ylabel('relative TEC (TECU)')
    title('Non-corrected non-calibrated slant TEC')
    savefig('_TEC1.png')
    close()

    #plot the peaks
    plot(param_sec, if400)
    plot(param_sec, if400_a)
    xlabel('time (s)')
    ylabel('FFT index (out of 16384)')
    title('400MHz peak indices')
    savefig('_peaks400.png')
    close()

    plot(param_sec, if150)
    xlabel('time (s)')
    ylabel('FFT index (out of 16384)')
    title('150MHz peak indices')
    savefig('_peaks150.png')
    close()

    #plot signal power
    sigpow150 = 10.0 * log10(sigpow150)
    plot(param_sec, sigpow150)
    xlabel('time (s)')
    ylabel('signal power (dB)')
    title('150MHz signal power')
    savefig('_sigpow150.png')
    close()

    sigpow400 = 10.0 * log10(sigpow400)
    plot(param_sec, sigpow400)
    xlabel('time (s)')
    ylabel('signal power (dB)')
    title('400MHz signal power')
    savefig('_sigpow400.png')
    close()

#########################################################################################

    pTEC = pTEC - correction * phase_sec #correct for the phase incursion
    #plot the corrected TEC
    plot(phase_sec,pTEC/1.0e16)
    xlabel('time (s)')
    ylabel('relative TEC (TECU)')
    title('Corrected non-calibrated slant TEC')
    savefig('_TEC2.png')
    close()

    (calibration,b) = polyfit(sec,pTEC,1)
    print 'calibration =', calibration/1.0e16, 'TECU'
    #plot the TEC vs. sec
    plot(sec,pTEC/1.0e16,'.')
    plot(sec,(sec*calibration+b)/1.0e16)
    xlabel('secant (1/rad)')
    ylabel('relative TEC (TECU)')
    title('Finding the calibration')
    savefig('_cali.png')
    close()

    pTEC = pTEC+calibration*sec[arg]-average(pTEC[arg-span:arg+span]) #calibrate for flat Vertical TEC
    #plot the calibrated TEC
    plot(phase_sec,pTEC/1.0e16)
    xlabel('time (s)')
    ylabel('relative TEC (TECU)')
    title('Calibrated slant TEC')
    savefig('_TEC3.png')
    close()

    #compare the secant and slant TEC
    plot(phase_sec,(pTEC-min(pTEC))/max(pTEC-min(pTEC)))
    plot(phase_sec,(sec-min(sec))/max(sec-min(sec)))
    xlabel('time (s)')
    ylabel('normalized values')
    title('matching the secant with the slant TEC')
    savefig('_comp.png')
    close()

    vTEC = pTEC*cos(ksi)
    print 'average TEC =', average(vTEC)/1.0e16, 'TECU'

    plot(phase_sec,vTEC/1.0e16)
    xlabel('time (s)')
    ylabel('relative TEC (TECU)')
    title('vertical TEC')
    savefig('_TEC4.png')
    close()

    print 'theoretical TEC minimum is at', phase_sec[ksi.argmin()], 's'
    print '     actual TEC minimum is at', phase_sec[pTEC.argmin()], 's'

    f_out2 = open('_peaks.txt', 'w')
    f_out2.write('time, 400MHz peak, 150MHz peak (Total %d lines)\n' % nparam)
    for ip in range(0, nparam):
        write_buf = '%g %g %g\n' % (param_sec[ip], if400[ip], if150[ip])
        f_out2.write(write_buf)
    f_out2.close

    f150.close
    f400.close

def _main():
    TECanalysis()
if __name__ == '__main__':
    _main()
