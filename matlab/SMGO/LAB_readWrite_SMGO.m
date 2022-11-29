function data = LAB_readWrite_SMGO(params,addr,duration,freq, lock_addr, jitter_addr,stepsize,t,kphi,dir,mean_)
    D = length(params);
    params_old = zeros(1,D);
    for i = 1:D
        data_str = doocsread(addr{i});
        params_old(i)=round(data_str.data,6);
    end
    counter=1:D;
    counter=counter(abs(params_old-params) > 10e-6);
    %disp(counter)
    for i=counter
       writeData(addr{i},params(i),params_old(i),stepsize(i),t)
       if i < length(counter)
            pause(1)
       end
    end
    pause(0.5)
    data = readData(addr{end},duration,freq,lock_addr,jitter_addr,kphi,dir,mean_);
end

function writeData(addr,param,param_old,stepsize,t)
    num_step=fix((param-param_old)/stepsize);
    disp(num_step)
    for i = 1:abs(num_step)
        if num_step > 0
            m=doocswrite(addr,param_old+i*stepsize);
        else
            m=doocswrite(addr,param_old-i*stepsize);
        end
        pause(t)
    end
    m=doocswrite(addr,param);
end

function [data]=readData(addr,duration,freq,lock_addr,jitter_addr,kphi,dir,mean_)
    times = 1/freq;
    time = duration;
    iter = ceil(time/times);
    jitter = zeros(iter,1);
    jitters = zeros(iter,length(jitter_addr));
    for i = 1:iter
         if ~checkLockStatus(lock_addr)
             count_str=load(dir+'/counter_unlock.mat');
             counter = count_str.counter;
             counter = counter + 1;
             save(dir+'/counter_unlock.mat','counter');
             st=tic;
            while ~checkLockStatus(lock_addr)
                pause(3)
                if toc(st) > 30
                    s = input("Can't lock system, set value manually by writing 'user' or continue by pressing 'enter': ",'s');
                    if strcmp(s,'user')
                        data = 20;
                        jit_str = load(dir+'/jitter_data','jitter_data');
                        jitter_data = jit_str.jitter_data;
                        jitter_data(end+1,:)=100*ones(1,length(jitter_addr)+1);
                        save(dir+'/jitter_data.mat','jitter_data')
                        return
                    else
                        st = tic;
                        continue;
                    end
                end
            end
         end
    
       data_str = doocsread(addr);
       jitter(i) = data_str.data*kphi;
       for l = 1:length(jitter_addr)
           temp_str = doocsread(jitter_addr{l});
           jitters(i,l) = temp_str.data;
       end
       pause(times)
    end
    if mean_
        data = mean(jitter);
    else
        data = median(jitter);
    end
    
    std_dev = std(jitter);
    jit_str = load(dir+'/jitter_data','jitter_data');
    jitter_data = jit_str.jitter_data;
    jitter_data(end+1,:) = [mean(jitters,1),std_dev];
    save(dir+'/jitter_data.mat','jitter_data');
    disp("measured value: " + num2str(data))
end

function b=checkLockStatus(lock_addr)
    c = 0;
    for i =1:length(lock_addr)
        c_str = doocsread(lock_addr{i});
        c = c + c_str.data;
    end
    if c == 0
        b = true;
    else
        b = false;
    end
end  

