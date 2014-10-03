'''
Created on Dec 22, 2013

@author: Sergio Hinojosa Rojas
'''

import numpy as np
from scipy.special import erfc
from scipy import interpolate, signal
import pywt
import glob
import os
import scipy
import shutil

class Sorter():
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.set_parameters()
        
    def set_parameters(self,
                w_pre = 20,                     # number of pre-event data points stored
                w_post = 44,                    # number of post-event data points stored
                detection = 'pos',              # type of threshold (pos, neg, both)
                stdmin = 5.00,                  # minimum threshold
                stdmax = 50,                    # maximum threshold
                interpolation = 'y',            # interpolation for alignment
                int_factor = 2,                 # interpolation factor
                detect_fmin = 300,              # high pass filter for detection (default 300)
                detect_fmax = 3000,             # low pass filter for detection (default 3000)
                sort_fmin = 300,                # high pass filter for sorting (default 300)
                sort_fmax = 3000,               # low pass filter for sorting (default 3000)
                max_spk = 20000,                # max. # of spikes before starting templ. match.
                template_type = 'center',       # nn, center, ml, mahal
                template_sdnum = 3,             # max radius of cluster in std devs.
                permut = 'y',                   # for selection of random 'par.max_spk' spikes before starting templ. match. 
                features = 'wav',               # choice of spike features (wav)
                inputs = 10,                    # number of inputs to the clustering (def. 10)
                scales = 4,                     # scales for wavelet decomposition
                mintemp = 0,                    # minimum temperature
                maxtemp = 0.301,                # maximum temperature (0.201)
                tempstep = 0.01,                # temperature step
                stab = 0.8,                     # stability condition for selecting the temperature
                SWCycles = 100,                 # number of montecarlo iterations (100)
                KNearNeighb = 11,               # number of nearest neighbors
                randomseed = 0,                 # if 0, random seed is taken as the clock value
                fname_in = 'tmp_data',          # temporary filename used as input for SPC
                fname = 'data_testfile',
                min_clus_abs = 20,              # minimum cluster size (absolute value)
                min_clus_rel = 0.005,           # minimum cluster size (relative to the total nr. of spikes)
                temp_plot = 'log',              # temperature plot in log scale (lin, log)
                force_auto = 'y',               # automatically force membership if temp>3.
                max_spikes = 5000,              # maximum number of spikes to plot.
                sr = 14400,                      # sampling frequency, in Hz.
                cluster_linux_dir = '/home/sergio/iibm/workspace2/NeuroDB/src/NeuroDB/signalProcessor'
                ):
        """ w_pre                 number of pre-event data points stored (default = 20)
            w_post                number of post-event data points stored (default = 44)
            detection             type of threshold (pos, neg, both) (default = 'pos')
            stdmin                minimum threshold (default = 5.00)
            stdmax                maximum threshold (default = 50)
            interpolation         interpolation for alignment (default = 'y')
            int_factor            interpolation factor (default = 2)
            detect_fmin           high pass filter for detection (default = 300)
            detect_fmax           low pass filter for detection (default = 3000)
            sort_fmin             high pass filter for sorting (default = 300)
            sort_fmax             low pass filter for sorting (default = 3000)
            max_spk               max.  of spikes before starting templ. match. (default = )
            template_type         nn, center, ml, mahal (default = 'center')
            template_sdnum        max radius of cluster in std devs. (default = 3 )
            permut                for selection of random 'par.max_spk' spikes before starting templ. match.  (default = 'y')
            features              choice of spike features (wav) (default = 'wav')
            inputs                number of inputs to the clustering (default = 10)
            scales                scales for wavelet decomposition (default = 4)
            mint                  minimum temperature (default = 0)
            maxtemp               maximum temperature (default = 0.301)
            tempstep              temperature step (default = 0.01)
            stab                  stability condition for selecting the temperature (default = 0.8)
            SWCycles              number of montecarlo iterations (100) (default = 100)
            KNearNeighb           number of nearest neighbors (default = 11)
            randomseed            if 0, random seed is taken as the clock value (default = 0)
            fname_in              temporary filename used as input for SPC (default = 'tmp_data')
            min_clus_abs          minimum cluster size (absolute value) (default = 20)
            min_clus_rel          minimum cluster size (relative to the total nr. of spikes) (default = 0.005)
            temp_plot             temperature plot in log scale (lin, log) (default = 'log')
            force_auto            automatically force membership if temp>3. (default = 'y')
            max_spikes            maximum number of spikes to plot. (default = 5000)
            sr                    sampling frequency, in Hz. (default = 14400) """
        
        self.w_pre =  w_pre
        self.w_post = w_post 
        self.detection =  detection
        self.stdmin = stdmin
        self.stdmax = stdmax
        self.interpolation = interpolation
        self.int_factor = int_factor
        self.detect_fmin = detect_fmin
        self.detect_fmax = detect_fmax
        self.sort_fmin = sort_fmin
        self.sort_fmax = sort_fmax
        self.max_spk = max_spk
        self.template_type = template_type
        self.template_sdnum = template_sdnum
        self.permut = permut
        self.features = features
        self.inputs = inputs
        self.scales = scales
        self.mintemp = mintemp
        self.maxtemp = maxtemp
        self.tempstep = tempstep
        self.stab = stab
        self.SWCycles = SWCycles
        self.KNearNeighb = KNearNeighb
        self.randomseed = randomseed
        self.fname_in = fname_in
        self.fname = fname
        self.min_clus_abs = min_clus_abs
        self.min_clus_rel = min_clus_rel
        self.temp_plot = temp_plot
        self.force_auto = force_auto
        self.max_spikes = max_spikes
        self.sr = sr
        self.cluster_linux_dir = cluster_linux_dir
    
    def __save_to_file(self, fname, data):        
        f = open (fname, "w")
        for row in data:
            for i in row:
                f.write(str(i))
                f.write(' ')
            f.write('\n')
            
        f.close()

    def __test_ks(self, x):
        x = x[np.invert(np.isnan(x))]
        n = len(x)
        x = np.sort(x)
        # Get cumulative sums
        yCDF = np.array(np.arange(n)+1) / float(n)
        # Remove duplicates; only need final one with total count
        notdup = np.append(np.diff(x), 1) > 0
        x_expcdf = x[notdup]
        y_expcdf = np.append(0, yCDF[notdup])
        
        # The theoretical CDF (theocdf) is assumed to be normal  
        # with unknown mean and sigma
    
        zScores  =  (x_expcdf - np.mean(x))/np.std(x, ddof=1)
        
        mu = 0
        sigma = 1
        theocdf = 0.5 * erfc(-(zScores-mu)/(np.sqrt(2)*sigma))
        
        delta1    =  y_expcdf[0:len(y_expcdf)-1] - theocdf   # Vertical difference at jumps approaching from the LEFT.
        delta2    =  y_expcdf[1:len(y_expcdf)] - theocdf   # Vertical difference at jumps approaching from the RIGHT.
        deltacdf  =  np.abs(np.concatenate((delta1, delta2), axis=0))
    
        return  np.max(deltacdf)
    
    def __wavedec(self, x , n, type):
        Lo_D = [0.7071,0.7071]
        Hi_D = [-0.7071,0.7071]
        c = []
        l = np.zeros([n+2])
        
        if (x == []):
            return
        
        l[len(l)-1] = len(x)
        for k in range(n):
            w = pywt.dwt(x,type)        # decomposition
            x = w[0]
            c =  np.concatenate((w[1],c),axis=0)            # store detail
            l[n-k] = len(w[1])    # store length
            
        c = np.concatenate((w[0],c),axis=0)
        l[0] = len(w[0])
    
        return c , l
        
    def __wave_features(self, spikes):
        nspk = len(spikes)
        ls = len(spikes[0])
        cc = np.zeros([nspk,ls])
        for i in range(nspk):
            c, l = self.__wavedec(spikes[i,:],4,'db1')
            cc[i,0:ls] = c[0:ls]
        
        sd = np.array([])
        for i in range(ls):             # KS test for coefficient selection   
            thr_dist = np.std(cc[:,i], ddof=1) * 3
            thr_dist_min = np.mean(cc[:,i]) - thr_dist
            thr_dist_max = np.mean(cc[:,i]) + thr_dist
            aux = cc[ [j for (j, val) in enumerate(cc[:,i]) if (val > thr_dist_min and val < thr_dist_max)],i ]
            if len(aux) > 10:
                ksstat = self.__test_ks(aux)
                sd = np.append(sd, ksstat)
            else:
                sd = np.append(sd, 0)
    
        ind = sd.ravel().argsort()
        coeff = ind[ls-1:ls-self.inputs-1:-1]
        
        inspk = np.zeros([nspk,self.inputs])
        for i in range(nspk):
            for j in range(self.inputs):
                inspk[i,j]=cc[i,coeff[j]]
        
        return inspk

    def __file2matrix(self, fid):
        f = fid.readlines()
        m =[]
        for row in f:
            row = row.strip()
            row = row.split(' ')
            row = [float(i) for i in row if i!=""]
            m.append(row)
            
        m = np.array(m)
        fid.close()
        return m

    def __run_cluster(self, fname_in, fname):
        dim = self.inputs
        
        # DELETE PREVIOUS FILES
        tmpfiles = glob.glob("*.dg_01*")
        if tmpfiles != []:
            for f in tmpfiles:
                os.remove(f)
        
        dat = open(fname_in)
        n = len(dat.readlines())
        dat.close()
        fid = open(fname + '.run','w')
        fid.write('NumberOfPoints: %s\n'%str(n))
        fid.write('DataFile: %s\n'%fname_in)
        fid.write('OutFile: %s\n'%fname)
        fid.write('Dimensions: %s\n'%str(dim))
        fid.write('MinTemp: %s\n'%str(self.mintemp))
        fid.write('MaxTemp: %s\n'%str(self.maxtemp))
        fid.write('TempStep: %s\n'%str(self.tempstep))
        fid.write('SWCycles: %s\n'%str(self.SWCycles))
        fid.write('KNearestNeighbours: %s\n'%str(self.KNearNeighb))
        fid.write('MSTree|\n')
        fid.write('DirectedGrowth|\n')
        fid.write('SaveSuscept|\n')
        fid.write('WriteLables|\n')
        fid.write('WriteCorFile~\n')
        if self.randomseed != 0:
            fid.write('ForceRandomSeed: %s\n'%str(self.randomseed))   
        fid.close()
        
        shutil.copy2(self.cluster_linux_dir+'/cluster_linux.exe', '/tmp/cluster_linux.exe')
        os.chmod('/tmp/cluster_linux.exe',664)
        os.system('chmod a+x /tmp/cluster_linux.exe')
        shutil.copy('%s.run'%fname, '/tmp/%s.run'%fname)
        #run_linux = './OpenSpikeSorter/cluster_linux.exe %s.run'%fname
        run_linux = '/tmp/cluster_linux.exe %s.run'%fname
        os.system(run_linux)
        
        fclu = open(fname + '.dg_01.lab')
        ftree = open(fname + '.dg_01')
        
        clu = self.__file2matrix(fclu)
        tree = self.__file2matrix(ftree)
        
        os.remove('%s.run'%fname)
        
        tmpfiles = glob.glob("*.mag") + glob.glob("*.edges") + glob.glob("*.param")
        
        if tmpfiles != []:
            for f in tmpfiles:
                os.remove(f)
                
        os.remove(fname_in)
        
        return clu, tree
    
    def __find_temp(self, tree,
                   num_temp = 30,
                   min_clus = 20):
        # Selects the temperature.
        
        aux  = np.diff(tree[:,4]);   # Changes in the first cluster size
        aux1 = np.diff(tree[:,5]);   # Changes in the second cluster size
        aux2 = np.diff(tree[:,6]);   # Changes in the third cluster size
        aux3 = np.diff(tree[:,7]);   # Changes in the third cluster size
        
        temp = 0;         # Initial value
        
        for t in range(num_temp):
            # Looks for changes in the cluster size of any cluster larger than min_clus.
            if ( aux[t] > min_clus or aux1[t] > min_clus or aux2[t] > min_clus or aux3[t] > min_clus ):    
                temp = t+1
                
        # In case the second cluster is too small, then raise the temperature a little bit 
        if (temp == 0 and tree[temp,5] < min_clus):
            temp = 1
        
        return temp  
        
    def do_clustering(self, spikes, index):
        num_temp = np.floor((self.maxtemp - self.mintemp)/self.tempstep)
        nspk = len(spikes)
        inspk = self.__wave_features(spikes)
        naux = np.minimum(self.max_spk,len(spikes))
        min_clus = np.maximum(self.min_clus_abs,self.min_clus_rel*naux)
    
        if self.permut == 'y':
            # GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
            if len(spikes[1,:]) > self.max_spk:
                # random selection of spikes for SPC 
                ipermut = np.random.permutation(range((len(inspk))))
                ipermut[naux+1:len(ipermut)] = []
                inspk_aux = inspk[ipermut,:]
            else:
                ipermut = np.random.permutation(range(len(inspk)))
                inspk_aux = inspk[ipermut,:];
            
    
            # INTERACTION WITH SPC
            self.__save_to_file(self.fname_in, inspk_aux)
            #save(handles.par.fname_in,'inspk_aux','-ascii');
            clu, tree = self.__run_cluster(self.fname_in, self.fname)
            temp = self.__find_temp(tree)
    
            # DEFINE CLUSTERS
            class1 = ipermut[np.argwhere(clu[temp,2:len(clu[0,:])]==0)]
            class1 = np.resize(class1,(1,len(class1)))[0]
            class2 = ipermut[np.argwhere(clu[temp,2:len(clu[0,:])]==1)]
            class2 = np.resize(class2,(1,len(class2)))[0]
            class3 = ipermut[np.argwhere(clu[temp,2:len(clu[0,:])]==2)]
            class3 = np.resize(class3,(1,len(class3)))[0]
            class4 = ipermut[np.argwhere(clu[temp,2:len(clu[0,:])]==3)]
            class4 = np.resize(class4,(1,len(class4)))[0]
            class5 = ipermut[np.argwhere(clu[temp,2:len(clu[0,:])]==4)]
            class5 = np.resize(class5,(1,len(class5)))[0]
            class0 = np.setdiff1d(range(len(spikes)), np.sort(np.concatenate((class1, class2, class3, class4, class5))))
    
            #whos class*
            
            # IF TEMPLATE MATCHING WAS DONE, THEN FORCE
            if (len(spikes)> self.max_spk or self.force_auto == 'y'):
                classes = np.zeros(len(spikes))
                if len(class1)>=min_clus:
                    classes[class1] = 1
                if len(class2)>=min_clus:
                    classes[class2] = 2
                if len(class3)>=min_clus:
                    classes[class3] = 3
                if len(class4)>=min_clus: 
                    classes[class4] = 4
                if len(class5)>=min_clus:
                    classes[class5] = 5
                    
                f_in  = spikes[classes!=0,:]
                f_out = spikes[classes==0,:]
                class_in = classes[np.argwhere(classes!=0),:]
                class_in = np.resize(class_in, (1,len(class_in)))[0]
            #class_out = force_membership_wc(f_in, class_in, f_out, handles)
            #classes(classes==0) = class_out
            #class0=find(classes==0)
            #class1=find(classes==1)
            #class2=find(classes==2)
            #class3=find(classes==3)
            #class4=find(classes==4)
            #class5=find(classes==5)
            cluster = np.zeros([nspk,2])
            #cluster[:,1]= np.resize(index,(len(index),1))
            cluster[:,1]=index
            num_clus = len(np.argwhere([len(class0), len(class1), len(class2), len(class3), len(class4), len(class5)] > min_clus))
            
            if (len(class0) > min_clus):
                max_spikes=np.minimum(len(class0),self.max_spikes)
                
            if len(class1) > min_clus:
                cluster[class1[:],0]=1
                
            if len(class2) > min_clus:
                cluster[class2[:],0]=2
                
            if len(class3) > min_clus:
                cluster[class3[:],0]=3
                
            if len(class4) > min_clus:
                cluster[class4[:],0]=4
                
            if len(class5) > min_clus:
                cluster[class5[:],0]=5
                
            #cluster_class = cluster
            np.savetxt('times_cluster.txt',cluster)
            
            if self.permut == 'n':
                np.savetxt('times_cluster.txt', cluster)
                np.savetxt('inspk.txt', inspk)
            else:
                #fid.write(par)
                np.savetxt('times_cluster.txt', cluster)
                np.savetxt('inspk.txt', inspk)
                np.savetxt('ipermut.txt', ipermut)
                
            return cluster, inspk
        

class Detector():
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.set_parameters()
        
    def set_parameters(self,
                         w_pre = 20,
                         w_post = 44,
                         detection = 'pos',
                         stdmin = 5.00,
                         stdmax = 50,
                         fmin_detect = 300,
                         fmax_detect = 3000,
                         fmin_sort = 300,
                         fmax_sort = 3000,
                         interpolation = 'y',
                         int_factor = 2,
                         sr = 14400,
                         min_ref_per = 1.5,
                         cluster_linux_dir = '/home/sergio/iibm/workspace2/NeuroDB/src/NeuroDB/signalProcessor'
                       ):
        ''' w_pre           number of pre-event data points stored.
            w_post          number of post-event data points stored.
            detection       type of threshold ('neg','both','pos')
            stdmin          minimum threshold.
            stdmax          maximum threshold.
            fmin_detect     high pass filter for detection.
            fmax_detect     low pass filter for detection.
            fmin_sort       high pass filter for sorting.
            fmax_sort       low pass filter for sorting.
            interpolation   interpolation for alignment.
            int_factor      interpolation factor.
            sr              sampling frequency, in Hz.
            min_ref_per     detector dead time. 
            ref             number of counts corresponding to the dead time'''
        
        self.w_pre = w_pre
        self.w_post = w_post
        self.detection = detection
        self.stdmin = stdmin
        self.stdmax = stdmax
        self.fmin_detect = fmin_detect
        self.fmax_detect = fmax_detect
        self.fmin_sort = fmin_sort
        self.fmax_sort = fmax_sort
        self.interpolation = interpolation
        self.int_factor = int_factor
        self.sr = sr
        self.min_ref_per = min_ref_per
        self.ref = np.floor(min_ref_per*sr/1000.0)
        self.cluster_linux_dir = cluster_linux_dir
        
    def __int_spikes(self, spikes):
                    
        ls = self.w_pre + self.w_post
        nspk = len(spikes[:,1])
        
        s = range(len(spikes[1,:]))
        #ints = range(1.0/int_factor,len(spikes[1,:]),1.0/int_factor)
        ints = np.arange(1.0/self.int_factor,len(spikes[1,:])+1.0/self.int_factor,1.0/self.int_factor)
        
        intspikes = np.zeros([1,len(ints)])
        spikes1 = np.zeros([nspk,ls])
        
        try:
            for i in range(nspk):
                if spikes[i,:].any():
                    intspikes = interpolate.spline(s,spikes[i,:],ints)
                    iaux = intspikes[self.w_pre*self.int_factor-3:self.w_pre*self.int_factor+6].argmax(0)
                    iaux = iaux + self.w_pre*self.int_factor
                    spikes1[i,0:self.w_pre] = intspikes[iaux-self.w_pre*self.int_factor-1:iaux-1:self.int_factor]
                    #spikes1[i,w_pre-1:-1:-1] = intspikes[iaux-int_factor:iaux-w_pre*int_factor-int_factor:-int_factor]
                    spikes1[i,self.w_pre:ls+1] = intspikes[iaux-1:iaux+self.w_post*self.int_factor-1:self.int_factor]
                else:
                    spikes1[i,:] = spikes[i,0:ls]
        except:
            print i
        return spikes1

    def __amp_detect(self, x):
        
        ref = np.floor(self.min_ref_per*self.sr/1000.0)
        
        # HIGH-PASS FILTER OF THE DATA
        (b,a) = signal.ellip(2, 0.1, 40, [self.fmin_detect*2.0/self.sr,self.fmax_detect*2.0/self.sr], btype='bandpass', analog=0, output='ba')
        xf_detect = signal.filtfilt(b, a, x)
        (b,a) = signal.ellip(2, 0.1, 40, [self.fmin_sort*2.0/self.sr,self.fmax_sort*2.0/self.sr], btype='bandpass', analog=0, output='ba')
        xf = signal.filtfilt(b, a, x)
        
        
        noise_std_detect = scipy.median(np.abs(xf_detect))/0.6745;
        noise_std_sorted = scipy.median(np.abs(xf))/0.6745;
       
        thr = self.stdmin * noise_std_detect        #thr for detection is based on detected settings.
        thrmax = self.stdmax * noise_std_sorted     #thrmax for artifact removal is based on sorted settings.
        
        # LOCATE SPIKE TIMES
        nspk = 0;
        xaux = np.argwhere(xf_detect[self.w_pre+1:len(xf_detect)-self.w_post-1-1] > thr) + self.w_pre + 1
        xaux = np.resize(xaux,len(xaux))
        xaux0 = 0;
        index = []
        for i in range(len(xaux)):
            if xaux[i] >= (xaux0 + ref):
            # after find a peak it begin search after ref over the last xaux
                iaux = xf[xaux[i]:xaux[i]+np.floor(ref/2.0)].argmax(0)    # introduces alignment
                nspk = nspk + 1
                index.append(iaux + xaux[i])
                xaux0 = index[nspk-1];
        
        # SPIKE STORING (with or without interpolation)
        ls = self.w_pre + self.w_post
        spikes = np.zeros([nspk,ls+4])
        xf = np.concatenate((xf,np.zeros(self.w_post)),axis=0)
        
        for i in range(nspk):                          # Eliminates artifacts
            if np.max( np.abs( xf[index[i]-self.w_pre:index[i]+self.w_post] )) < thrmax :
                spikes[i,:] = xf[index[i]-self.w_pre-1:index[i]+self.w_post+3]
     
        aux = np.argwhere(spikes[:,self.w_pre] == 0)       #erases indexes that were artifacts
        if len(aux) != 0:
            aux = aux[0][0]
            spikes = np.delete(spikes, aux, axis = 0)
            index = np.delete(index,aux)
 
        if self.interpolation == 'y':
            # Does interpolation
            spikes = self.__int_spikes(spikes)

        return spikes, thr, index
        
    def get_spikes(self, analogSignal):
        spikes, thr, index = self.__amp_detect(analogSignal)
        
        return spikes, index, thr

if __name__ == '__main__':
    pass