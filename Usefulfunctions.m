classdef Usefulfunctions
    methods (Static) % General functions
        function b=iseven(x)
            b=Usefulfunctions.isint(x/2);
        end
        function b=isodd(x)
            b=Usefulfunctions.isint(x) && ~Usefulfunctions.isint(x/2);
        end
        function b=isint(x)
            b= (length(x)==1) && (isreal(x)) && (mod(x,1)==0) ;
        end
        function addtexttofile(text,file)
            
            fid=fopen(file,'a+');
            fprintf(fid,'%s\n',text);
            fclose(fid);
        end
        function dispnum(varargin)
            
            s=[];
            for i=1:nargin
                s=[s num2str(varargin{i}) ' ']; %#ok<*AGROW>
            end
            disp(s);
        end
        function disptextnum(varargin)
            
            s=[];
            for i=1:2:nargin-1
                s=[s varargin{i} num2str(varargin{i+1})];
            end
            if nargin/2~=round(nargin/2)
                s=[s varargin{nargin}];
            end
            disp(s);
        end
        function dirdate
            m=dir;
            a=zeros(length(m),1);
            for i=1:length(m)
                a(i)=m(i).datenum;
            end
            [~,ind]=sort(a);
            disp('------Directories------');
            for i=1:length(a)
                j=ind(i);
                if m(j).isdir
                    disp(m(j).name);
                end
            end
            disp('------Files------');
            for i=1:length(a)
                j=ind(i);
                if ~m(j).isdir
                    disp([m(j).name, char(32*ones(1,30-length(m(j).name))),m(j).date]);
                end
            end
        end
    end
    methods (Static) %Signal Processing fuctions
        function resetrandom
            
            s = RandStream('mt19937ar','Seed',0);
            RandStream.setGlobalStream(s);
        end
        function E=energy(x)
            E=sum(Usefulfunctions.abs2(x));
        end
        function y=abs2(x)
            y=abs(x).^2;
        end
        function x=randn_c( rows, cols, thresh )
            %RANDN_C Complex normal random numbers, neg
            if nargin<3, thresh=0; end
            if nargin<2
                if nargin<1
                    rows=1; cols=1;
                else
                    if length(rows)==2
                        cols=rows(2);
                        rows=rows(1);
                    else
                        error('bad arguments');
                    end
                end
            end
            
            if thresh<=0
                x=sqrt(thresh^2-2*log(rand(rows,cols))).*exp(1i*2*pi*rand(rows,cols));
            else
                x=sqrt(-2*log(1-(1-exp(-thresh^2/2))*rand(rows,cols))).*exp(1i*2*pi*rand(rows,cols));
            end
            
        end
        function y=delay( x, n , OS)
            %DELAY(X,N) Delays the vector x by n steps. If n is ommitted it is set to
            % 1.
            % Updated 2017-07-06: bug found. Changed from t+frac to to t-frac
            % Updated 2018-12-23: allowing for vector of delays n
            % Updated 2020-05-12: allowing for matrix x
            % (c) Thomas Eriksson 2017-2020
            
            
            if nargin<3, OS=1; end
            if nargin<2, n=1; end
            
            cols=size(x,2);
            
            y=zeros(size(x,1),length(n)*cols);
            for ind=1:length(n)
                nint=round(n(ind));
                y=circshift(x,nint);
                if n(ind)>=0
                    y(1:nint,:)=zeros(nint,cols);
                else
                    y(end+nint+1:end,:)=zeros(abs(nint),cols);
                end
                frac=n(ind)-nint;
                if frac~=0 % check if integer
                    if OS>=2  % simple version if there is some oversampling
                        O=5;
                        h=rc((-O:O)'-frac,0.5);
                    else
                        O=512;
                        h=sinc((-O:O)'-frac).*hanning(2*O+1);
                    end
                    y=[y;zeros(O,cols)];
                    y=filter(h,1,y);
                    y=y(end-length(x)+1:end,:);
                end
                y(:,(ind-1)*cols+1:ind*cols)=y;
            end
        end
        function [y,p]=limiter(x,limit, limithi)
            % y=limiter(x, [limit], [limithi] )
            % Limits abs(x) to limit (defailt limit=1)
            if nargin<2
                limit=1;
            end
            
            y=x;
            if nargin<3
                i=find(abs(x(:))>limit);
                y(i)=limit*x(i)./abs(x(i));
            else
                i= x(:)<limit;
                y(i)=limit;
                i= x(:)>limithi;
                y(i)=limithi;
            end
            p=length(i)/length(x);
        end
        function x=rect(t)
            
            x=zeros(size(t));
            x(abs(t)<0.5)=1;
        end
        function x=tri(N)
            
            t=linspace(-1,1,N+2);
            x=triangle(t(2:end-1));
        end
        function x=impulse(L,pos)
            
            if nargin<2
                pos=1;
            end
            x=zeros(L,1);x(pos)=1;
            
        end
        function x=kroneckerdelta(t)
            
            x=zeros(size(t));
            x(t==0)=1;
        end
        function [out,alpha]=snr(x,y,type)
            % SNR(X,Y,type)
            % If type is 'rescaled' then y magnitude is chosen optimal.
            
            if nargin<3
                type='1';
            end
            
            switch type
                case 'rescaled'
                    alpha=(y(:)'*x(:))/sum(abs(y(:)).^2);
                case '1'
                    alpha=1;
                otherwise
                    error('wrong last argument');
            end
            out=10*log10(var(x(:))/var(x(:)-alpha*y(:)));
        end
        function s=spec(inputsignal,N,fs)
            % SPEC Double-sided spectrum, normalized frequency
            % SPEC(X) plots the spectrum of the vector X.
            % SPEC(X,N) plots the spectrum of X with N points resolution.
            % SPEC([X,Y,Z],N) plots the spectrum of X, Y and Z in the same plot.
            
            if nargin<3, fs=1; end
            if nargin<2, N=1024; end
            
            if size(inputsignal,2)>size(inputsignal,1)
                inputsignal=inputsignal';
            end
            if (length(inputsignal)<N)
                inputsignal=[zeros(floor((N-length(inputsignal))/2),size(inputsignal,2));inputsignal;zeros(ceil((N-length(inputsignal))/2),size(inputsignal,2))];
            end
            spectrumplotcolor=[0 0 1 ; 1 0 0 ; 0 1 0 ; 0 0 0 ; 1 1 0 ; 1 0 1; 0 1 1; 0.6 0.6 0.6];
            spectrumplotcolor=[spectrumplotcolor;0.6*spectrumplotcolor];
            legendstrings=[];
            for k=1:size(inputsignal,2)
                x=inputsignal(:,k);
                N=min(N,length(x)-1);
                [X,w]=pwelch(x*sqrt(2*pi),hanning(N),[],'centered');
                if nargout==1
                    s(:,k)=10*log10(X);
                else
                    plot(w/2/pi*fs,10*log10(X),'color',spectrumplotcolor(k,:),'linewidth',2);
                    xlabel('Normalized frequency');
                    ylabel('Spectrum [dB]');
                    grid on
                    hold on
                end
                legendstrings=[legendstrings;num2str(k)];
            end
            if nargout~=1
                axis tight;
                legend(legendstrings);
                hold off
            end
        end
        function y=wrap(x,f)
            if nargin<2
                f=2*pi;
            end
            t=floor((x+f/2)/f);
            y=x-f*t;
        end
        function t=time(x)
            
            if Usefulfunctions.isint(x)
                t=(1:x)';
            else
                t=(1:length(x))';
            end
        end
        function [W,fc,linf,SNR]=bandwidth_estimate(x)
            
            X=spec(x);
            Xl=10.^(spec(x)/10);
            f=(1:length(X))';f=f-mean(f);
            p=zeros(4,1);
            p(2)=f'*Xl/sum(Xl);
            m1=sum(Xl);
            m2=sum(Xl.^2);
            p(3)=m1^2/m2;
            
            
            H=[rect((f-p(2))/p(3)) ones(size(f))];
            th=pinv(H)*X;
            p(1)=th(1);
            p(4)=th(2);
            
            W=fminsearch(@(W) sum(abs2(Xl-10.^( (p(1).*rect((f-p(2))./W)+p(4))/10  ))),p(3));
            popt=fminsearch(@(p) sum(abs2(X-p(1).*rect((f-p(2))./W)-p(4))),p);
            
            W=W/length(X);
            fc=popt(2)/length(X);
            linf=popt(4);
            SNR=popt(1);
            
            if SNR<3
                W=1;
                SNR=inf;
                linf=-inf;
            end
        end
        function y=dsm(x)
            
            q=0;
            esum=0;
            y=zeros(size(x));
            for n=1:length(x)
                e=x(n)-q;
                esum=esum+e;
                q=sign(esum+eps);
                y(n)=q;
            end
        end
        function GMM=em(X, NOMIXTURES, NOITER, INITMETHOD, STOPVALUE)
            % EM    The EM algorithm for optimizing Gaussian mixture models (GMMs) with diagonal covariances
            %   The function can be called in two ways:
            %   GMM = EM(X,NOMIXTURES) trains a GMM with NOMIXURES mixtures
            %       X is the training matrix, size d*s, where each column is a training
            %       vector.
            %   GMM = EM(X,GMM) trains a GMM with parameters as in the input GMM
            %       GMM should contain fields w (1*M), m (d*M), sigma2 (d*M), and LL (1*1).
            % Thomas Eriksson, 2003
            % Email: thomase@s2.chalmers.se
            
            if nargin<3
                NOITER=200;  % maximum number of iterations
            end
            if nargin<4
                INITMETHOD='vq';
            end
            if nargin<5
                STOPVALUE=1e-3;
            end
            if NOITER<=0
                NOITER=200;
            end
            
            
            d=size(X,1); % dimension
            T=size(X,2); % number of training vectors
            if isfield(NOMIXTURES,'w')
                w=NOMIXTURES.w;
                m=NOMIXTURES.m;
                sigma2=NOMIXTURES.sigma2;
                M=length(w);
            else
                M=NOMIXTURES; % number of mixtures
                
                %rand('state',0);
                % initialize the output variables
                w=ones(1,M)/M;
                switch INITMETHOD
                    case 'vq'
                        g=sqrt(var(X,2));
                        m=vqtrain(X./repmat(g,1,size(X,2)),M,100000,0.5);
                        m=m.*repmat(g,1,size(m,2));
                    case 'random'
                        m=X(:,ceil(rand(1,M)*T));
                    case 'kmeans'
                        g=sqrt(var(X'))';
                        [~,m]=kmeans((X./repmat(g,1,size(X,2)))',M,'maxiter',10,'start','sample');m=m';
                        m=m.*repmat(g,1,size(m,2));
                end
                sigma2=repmat(sum((X-repmat(mean(X,2),1,length(X))).^2,2)/length(X),1,M);
            end
            
            LL=0;
            tic;
            LLold=-9999;
            n=0;
            while n<NOITER
                p_aposteriori=zeros(T,M);
                for i=1:M
                    s=zeros(1,length(X));
                    for j=1:d
                        s=s-0.5/sigma2(j,i) *  (   X(j,:)-m(j,i) ).^2   ;
                    end
                    p_aposteriori(:,i)=(  w(i)/(2*pi)^(d/2)./sqrt(prod(sigma2(:,i)))*exp(s)  )';
                end
                temp=sum(p_aposteriori,2);
                p_aposteriori=p_aposteriori./repmat(temp,1,M);
                w=sum(p_aposteriori);
                m=(X*p_aposteriori)./repmat(w,d,1);
                sigma2=(X.^2*p_aposteriori)./repmat(w,d,1)-m.^2;
                
                maxes=max(sigma2');
                for i=1:size(sigma2,1)
                    val=0.001*maxes(i);
                    for j=1:size(sigma2,2)
                        if sigma2(i,j)<val
                            disp(['Normalizing variance in dim ' num2str(i) ', mixture ' num2str(j)]);
                            sigma2(i,j)=val;
                        end
                    end
                end
                w=w/T;
                LL=sum(log(temp))/T;
                
                disp([num2str(n) ': LL = ', num2str(LL) ' ' num2str(toc)]);
                if LL-LLold<STOPVALUE
                    break;
                end
                LLold=LL;
                n=n+1;
                tic;
            end
            GMM.w=w;
            GMM.m=m;
            GMM.sigma2=sigma2;
            GMM.LL=LL;
            
        end
        function GMM=emfull(X, NOMIXTURES, NOITER, INITMETHOD)
            % EM    The EM algorithm for optimizing Gaussian mixture models (GMMs) with full covariances
            %   GMM = EM(X,NOMIXTURES) trains a GMM with NOMIXURES mixtures
            %       X is the training matrix, size d*s, where each column is a training
            %       vector.
            %   GMM = EM(X,NOMIXTURES,NOITER,INITMETHOD) uses maximally NOITER
            %   iterations. INITNETHOD can be 'random','vq', or 'kmeans'.
            %
            %   The function can alse be called with a gmm as starting point:
            %   GMM = EM(X,GMM) trains a GMM with parameters as in the input GMM
            %       GMM should contain fields w (1*M), m (d*M), sigma2 (d*M), and LL (1*1).
            % Thomas Eriksson, 2004
            % Email: thomase@s2.chalmers.se
            
            if nargin<3
                NOITER=200;  % maximum number of iterations
            end
            if nargin<4
                INITMETHOD='vq';
            end
            if NOITER<=0
                NOITER=200;
            end
            
            
            d=size(X,1); % dimension
            T=size(X,2); % number of training vectors
            if isfield(NOMIXTURES,'w')
                w=NOMIXTURES.w;
                m=NOMIXTURES.m;
                sigma2=NOMIXTURES.sigma2;
                M=length(w);
            else
                M=NOMIXTURES; % number of mixtures
                
                %rand('state',0);
                % initialize the output variables
                w=ones(1,M)/M;
                switch INITMETHOD
                    case 'vq'
                        m=vqtrain(X,M,100000,0.5);
                    case 'random'
                        m=X(:,ceil(rand(1,M)*T));
                    case 'kmeans'
                        [~,m]=kmeans(X',M,'maxiter',10,'start','sample');m=m';
                end
                A=(X-repmat(mean(X,2),1,length(X)));
                A=A*A'/length(X);
                sigma2=zeros(size(A,1),size(A,2),M);
                for i=1:M
                    sigma2(:,:,i)=A;
                end
            end
            
            LL=0;LLold=-9999;
            n=0;
            tic;
            while n<NOITER
                p_aposteriori=zeros(T,M);
                for i=1:M
                    A=(X-repmat(m(:,i),1,T))';
                    p_aposteriori(:,i)=w(i)/(2*pi)^(d/2)/sqrt(det(sigma2(:,:,i)))*exp(  -0.5*( sum(A*inv(sigma2(:,:,i)).*A,2) ) );
                end
                temp=sum(p_aposteriori,2);
                p_aposteriori=p_aposteriori./repmat(temp,1,M);
                w=sum(p_aposteriori);
                m=(X*p_aposteriori)./repmat(w,d,1);
                
                for i=1:M
                    sigma2(:,:,i)=X*(X'.*repmat(p_aposteriori(:,i),1,d))/w(i)-m(:,i)*m(:,i)';
                    if cond(sigma2(:,:,i))>500
                        [V,D] = eig(sigma2(:,:,i));
                        D = diag(max(diag(D),max(diag(D))/500));
                        sigma2(:,:,i) =V*D*V';
                    end
                end
                w=w/T;
                LL=sum(log(temp))/T;
                
                disp([num2str(n) ': LL = ', num2str(LL) ' ' num2str(toc)]);
                if LL-LLold<1e-3
                    break;
                end
                LLold=LL;
                n=n+1;
                tic;
            end
            
            GMM.w=w;
            GMM.m=m;
            GMM.sigma2=sigma2;
            GMM.LL=LL;
            
            
            %%%%%%%%%%%%%%% Below is the old, slower algorithm (but easier to
            %%%%%%%%%%%%%%% understand :-)
            % while n<NOITER
            %     for i=1:M
            %         normweight(1,i) = w(i)/(2*pi)^(d/2)/sqrt(det(sigma2(:,:,i)));
            %         normexp(:,:,i)  = -0.5*inv(sigma2(:,:,i));
            %     end
            %     wsum=zeros(1,M);
            %     msum=zeros(d,M);
            %     sigma2sum=zeros(d,d,M);
            %     LL=0;
            %     for t=1:T
            %         x1=X(:,t);
            % %        x2=diag(x1.^2);
            %         x2=x1*x1';
            %         for i=1:M
            %             p_aposteriori(i)=normweight(i)*exp(  (x1-m(:,i))'*normexp(:,:,i)*(x1-m(:,i)));
            %         end
            %         fx=sum(p_aposteriori);
            %         p_aposteriori=p_aposteriori/fx;
            %         LL=LL+log(fx);
            %         wsum=wsum+p_aposteriori;
            %         msum=msum+x1*p_aposteriori;
            %         for i=1:M
            %             sigma2sum(:,:,i)=sigma2sum(:,:,i)+x2*p_aposteriori(i);
            %         end
            %     end
            %      for i=1:M
            %         m(:,i)=msum(:,i)/wsum(i);
            % %        sigma2(:,:,i)=sigma2sum(:,:,i)./wsum(i)-diag(m(:,i).^2);
            %         sigma2(:,:,i)=sigma2sum(:,:,i)./wsum(i)-m(:,i)*m(:,i)';
            %     end
            %     w=wsum/T;
            %     LL=LL/T;
            %     disp([num2str(n) ': LL = ', num2str(LL) ' ' num2str(toc)]);
            %     if LL-LLold<1e-3
            %         break;
            %     end
            %     LLold=LL;
            %     n=n+1;
            %     tic;
            % end
        end
        function y=fftbpfilter(x,fc,bw)
            
            IND=round(length(x)*fc+length(x)/2);
            BWI=round(length(x)*bw);
            Y=fftshift(fft(x));
            Yf=eps*ones(size(Y));
            %Yf(IND-BWI:IND+BWI)=Y(IND-BWI:IND+BWI).*hanning(2*BWI+1);
            Yf(IND-BWI:IND+BWI)=Y(IND-BWI:IND+BWI);
            y=ifft(ifftshift(Yf));
            
        end
        function y=fftbsfilter(x,fc,bw)
            
            IND=round(length(x)*fc+length(x)/2);
            BWI=round(length(x)*bw);
            Y=fftshift(fft(x));
            Yf=Y;
            f=linspace(-1,1,2*BWI+1);
            r=(5*triangle(f)-4*triangle(f*5/4))';
            
            Yf(IND-BWI:IND+BWI,:)=Yf(IND-BWI:IND+BWI,:).*(1-r+1e-4);
            y=ifft(ifftshift(Yf));
            
            
        end
        function y=Q(x)
            
            out = 0.5*erfc(x/sqrt(2));
        end
        function y=sincinterpolate(x,t,N)
            % Y = SINCINTERPOLATE(X,T,N)
            % Computes the signal at arbitrary time positions T.
            % X is assumed to have a sampling interval of 1, with the first entry in X at time 1
            % and the last at time length(X).
            % T is a vector of new desired sampling instants.
            % The elements of T should normally be an increasing sequence with values
            % in the interval 1 <= T(k) < T(k+1) <= length(x).
            % N is the window size. (default 64)
            % (c) Thomas Eriksson 2007
            
            if nargin<3, N=64; end  % The length of the interpolation is 2*N+1, see calc of p1 and p2 below.
            
            if size(x,1)<=1
                x=x.';
            end
            t=t(:);
            y=zeros(length(t),size(x,2));
            for k=1:length(t)
                p=round(t(k));
                p1=max(1,min(p-N,size(x,1)-1));
                p2=min(max(p+N,2),size(x,1));
                h=0.5*(1+cos(pi/(p2-p1)*((p1:p2)-t(k))));
                y(k,:)=(h.*sinc((p1:p2)-t(k)))*x(p1:p2,:);
            end
            
        end
        function y=sincinterpolate2(x,t)
            % Y = SINCINTERPOLATE2(X,T)
            % THIS FUNCTION IS INTENDED FOR USE WITH OVERSAMPLED X
            % IT USES RC(0.5) PULSES INSTEAD OF SINC
            % Computes the signal at arbitrary time positions T.
            % X is assumed to have a sampling interval of 1, with the first entry in X at time 1
            % and the last at time length(X).
            % T is a vector of new desired sampling instants.
            % The elements of T should normally be an increasing sequence with values
            % (c) Thomas Eriksson 2017
            
            N=8;
            if size(x,1)<=1
                x=x.';
            end
            t=t(:);
            
            x=[zeros(N,size(x,2));x;zeros(N,size(x,2))];
            t=t+N;
            y=zeros(length(t),size(x,2));
            p1=round(t-N);
            p2=round(t+N);
            
            for k=1:length(t)
                t2=(p1(k):p2(k))-t(k)+eps;
                %    y(k,:)=1/2/pi*(sin(1.5*pi*t2)+sin(0.5*pi*t2))./(t2.*(1-t2.^2))*x(p1(k):p2(k),:);
                %    y(k,:)=rc(t2,0.5)*x(p1(k):p2(k),:);
                y(k,:)=(sinc(0.75*t2).^2-1/9*sinc(t2/4).^2)*x(p1(k):p2(k),:);
            end
            
            
        end
        function y=sincresample(x,f)
            % Y = SINCRESAMPLE(X,F)
            % Resamples X by a factor of F.
            % (c) Thomas Eriksson 2008
            
            assert(size(x,1)>size(x,2),'sincresample operates on columns');
            
            %k=(size(x,1)-1)/(f*size(x,1)-1);
            t=1:1/f:(size(x,1)+(f-1)/f);
            y=sincinterpolate([x(1,:);x],t+1);
        end
    end
    methods (Static) % Amplifier/DPD functions
        function out=acpr(x)
            
            x=x(20:end-19,:);
            
            channelBW=bandwidth_estimate(x);
            channelPos=channelBW*0.6;
            
            
            X=fft(x);
            X=fftshift(X);
            SX=abs(X).^2;
            L=length(x);
            M=L/2+1;
            Pmain=sum(SX(round(M-channelBW*L/2):round(M+channelBW*L/2)));
            Pside1=sum(SX(round(M+channelPos*L):round(M+channelPos*L+channelBW*L)));
            Pside2=sum(SX(round(M-channelPos*L-channelBW*L):round(M-channelPos*L)));
            Pside=max(Pside1,Pside2);
            
            out=10*log10(Pside/Pmain);
            
            
            if 0
                SXdB=10*log10(SX);
                plot(SXdB)
                line([round(M-channelBW*L/2),round(M-channelBW*L/2)],[min(SXdB),max(SXdB)],'color','r','linewidth',3)
                line([round(M+channelBW*L/2),round(M+channelBW*L/2)],[min(SXdB),max(SXdB)],'color','r','linewidth',3)
                line([round(M-channelPos*L-channelBW*L),round(M-channelPos*L-channelBW*L)],[min(SXdB),max(SXdB)],'color','g','linewidth',3)
                line([round(M-channelPos*L),round(M-channelPos*L)],[min(SXdB),max(SXdB)],'color','g','linewidth',3)
                line([round(M+channelPos*L+channelBW*L),round(M+channelPos*L+channelBW*L)],[min(SXdB),max(SXdB)],'color','g','linewidth',3)
                line([round(M+channelPos*L),round(M+channelPos*L)],[min(SXdB),max(SXdB)],'color','g','linewidth',3)
            end
        end
        function amam(x,y,varargin)
            clf
            plot(abs(x),abs(y),'.');
            hold on
            for arg=1:2:nargin-2
                plot(abs(varargin{arg}),abs(varargin{arg+1}),'.');
            end
            ylabel('|y|');
            xlabel('|x|');
            
            title('AM/AM plot');
            hold off
        end
        function ampm(x,y)
            plot(abs(x),wrap(angle(y)-angle(x)),'.');
        end
        function [alpha,dist]=bussgang(x,y)
            alpha=mean(conj(y).*x)/mean(conj(x).*x);
            dist=y-alpha*x;
        end
        function out=cubicmetric(x)
            % cubic metric
            
            x=x./sqrt(var(x));
            out=(10*log10(var(x.^3))- 1.5237)/1.85;
        end
        function A=dpdid_ila(pa, model, x)
            % B.name='gmp'; B.P=7; B.M=1; B.G=1;
            % A=dpdid_ila(@pa, B, x)
            
            % find the gain
            y=pa(0.01*x);
            gain=sqrt(energy(y)/energy(0.01*x));
            y=pa(x/gain);
            gain=sqrt(energy(y)/energy(x/gain));
            
            
            % a first round
            u=x/gain;
            y=pa(u);
            A=identify_model(y,u,model);
            
            % a second round
            u=run_model(x,A);
            y=pa(u);
            A=identify_model(y,u,model);
        end
        function [y,u]=dpdid_ilc(ydesired, pa, varargin)
            % allowed extra parameters are 'gain','xmax','iterations', all
            % with an additional numerical argument.
            
            % find the gain
            gain=-1;
            xmax=inf;
            
            norounds=25;
            arg=1;
            while arg<=nargin-2
                switch varargin{arg}
                    case 'gain'
                        gain=varargin{arg+1};
                        arg=arg+1;
                    case 'xmax'
                        xmax=varargin{arg+1};
                        arg=arg+1;
                    case 'iterations'
                        norounds=varargin{arg+1};
                        arg=arg+1;
                    otherwise
                        disp(varargin{arg});
                        error('Nonexisting option');
                end
                arg=arg+1;
            end
            
            if (gain<=0)
                y=pa(0.001*ydesired);
                gain=abs(pinv(0.001*ydesired)*y)
            end
            
            
            
            gammastart=1/2/gain;
            %gamma=1/sqrt(energy(y)/energy(u));
            % a first round
            u=ydesired/gain;
            u=Usefulfunctions.limiter(u,xmax);
            y=pa(u);
            
            disp(['SNR: ' num2str(Usefulfunctions.snr(ydesired,y,'rescaled')) '    DBM: ' num2str(Usefulfunctions.getdbm(y))]);
            
            for rounds=1:norounds
                gamma=(1-(rounds-1)/norounds)*gammastart;
                %gamma=gammastart;
                e=ydesired-y;
                u=u+gamma*e;
                
                us=sort(abs(u));
                a=us(round(0.9999*length(us)));
                a=min(a,xmax);
%                u=Usefulfunctions.limiter(u,a);
                
                y=pa(u);
                disp(['SNR: ' num2str(Usefulfunctions.snr(ydesired,y)) '    DBM: ' num2str(Usefulfunctions.getdbm(y))]);
            end
        end
        function x2=fftresample(x,f)
            % Y = fftresample(X,F)
            % resamples X by a factor of F (up or down)
            % f can be a real value.
            % (c) Thomas Eriksson 2015
            
            assert(size(x,1)>size(x,2),'fftresample operates on columns');
            L1=length(x);
            N1=pow2(nextpow2(L1)+1);
            x=[x;zeros(N1-L1,size(x,2))];
            
            L=length(x);
            N=round(L*f);
            X=fft(x);
            % The cases: N>L, N=L, N<L, and N odd, N even
            if N<L
                if floor(N/2)==N/2
                    X2=[X(1:N/2);X(N/2+1)+X(end-N/2+1);X(end-N/2+2:end)];
                else
                    X2=[X(1:ceil(N/2));X(end-floor(N/2)+1:end)];
                end
            else
                if N>L
                    if floor(L/2)==L/2
                        X2=[X(1:L/2);X(L/2+1)/2;zeros(N-L-1,1);X(end-L/2+1)/2;X(end-L/2+2:end)];
                    else
                        X2=[X(1:ceil(L/2));zeros(N-L,1);X(end-floor(L/2)+1:end)];
                    end
                else
                    X2=X;
                end
            end
            x2=f*ifft(X2);
            
            
            x2=x2(1:round(L1*f));
        end
        function y=ghorbani(x)
            % [1] A. Ghorbani, and M. Sheikhan, "The effect of Solid State Power Amplifiers (SSPAs) Nonlinearities on MPSK and M-QAM Signal Transmission", Sixth Int'l Conference on Digital Processing of Signals in Comm., 1991, pp. 193-197.
            
            x1=0;
            x2=0;
            x3=0;
            x4=0;
            y1=0;
            y2=0;
            y3=0;
            y4=0;
            
            
            r=abs(x);
            f=angle(x);
            
            r2=x1*r^x2./(1+x3*r^x2)+x4*r;
            f2=f+y1*r^y2./(1+y3*r^y2)+y4*r;
            y=r2.*exp(1i*f2);
        end
        function y=poly5(x)
            
            alpha=10/9-sqrt(20/9*10^-0.05-80/81); %This factor gives 1dB compression at x=1
            beta=9*alpha^2/20; % gives minimum derivative=0, at x=sqrt(2/3/alpha)
            
            r=abs(x);
            f=angle(x);
            r2 = r.*(1-alpha*r.^2+beta*r.^4);
            
            y=r2.*exp(1i*f);
        end
        function y=rapp(x)
            % [1] C. Rapp, "Effects of HPA-Nonlinearity on a 4-DPSK/OFDM-Signal for a Digitial Sound Broadcasting System", in Proceedings of the Second European Conference on Satellite Communications, Liege, Belgium, Oct. 22-24, 1991, pp. 179-184.
            
            A0=1;
            p=3.0;
            v=1.0;
            A=abs(x);
            f=angle(x);
            
            A2=A*v./(1+(A*v/A0).^(2*p)).^(1/2/p);
            
            y=A2.*exp(1i*f);
        end
        function y=saleh(x)
            % [1] Saleh, A.A.M., "Frequency-independent and frequency-dependent nonlinear models of TWT amplifiers," IEEE Trans. Communications, vol. COM-29, pp.1715-1720, November 1981.
            
            alpha_am=2.0;
            beta_am=1.0;
            alpha_pm=pi/3;
            beta_pm=1.0;
            
            r=abs(x);
            f=angle(x);
            
            r2=alpha_am*r./(1+beta_am*r.^2);
            f2=f+alpha_pm*r.^2./(1+beta_pm*r.^2);
            y=r2.*exp(1i*f2);
        end
        function y=softlimiter(x,p)
            if nargin<2
                p=1;
            end
            y=(atan(abs(x).^(1/p))*2/pi).^p.*exp(1i*angle(x));
        end
        function out=metric(x)
            
            sigmax=sqrt(mean(abs(x).^2));
            
            out=(mean(abs(x).^6)/sigmax^6-(mean(abs(x).^4))^2/sigmax^8);
        end
        function out=papr(x,part,OS)
            
            if nargin<3 OS=1; end
            if nargin<2 part=1; end
            
            if OS>1
                x=resample(x(:),OS,1);
            end
            x=sort(abs(x).^2);
            maxval=x(round(part*length(x)));
            meanval=mean(x);
            
            out=10*log10(maxval/meanval);
            
        end
        function out=rawcubicmetric(x)
            % raw cubic metric
            
            x=abs(x);
            x=x./rms(x);
            out=20*log10(rms(x.^3));
        end
        function out=getdbm(x)
            % columnwise dBm
            out=10*log10(mean(abs(x).^2,1)/50/1e-3);
        end
        function out=setdbm(x,dnew)
            
            if Usefulfunctions.getdbm(x)>-100
                out=x.*repmat(10.^((dnew-Usefulfunctions.getdbm(x))/20),size(x,1),1);
            else
                out=x;
            end
        end
        function out=nmse(x,y)
            % nmse(x,y)
            % Remove samples from beginning and end, due to time alignment issues.
            alpha=(y(10:end-9)'*x(10:end-9))/sum(abs(y(10:end-9)).^2);
            out=10*log10(var(alpha*y(10:end-9)-x(10:end-9),1)/var(x(10:end-9),1));
        end
        function [y2,timeadjustment]=timealign(x,y,mode)
            % TIMEALIGN - function to align one signal to another
            %
            % y2=timealign(x,y) shifts the sampling instants of y such that x and y2
            % are aligned, to a subsample resolution.
            % y2=timealign(x,y,'phasetracking') also includes tracking and
            % compensation of phase noise in y.
            % (c) Thomas Eriksson 2015
            % 2017: Added phase noise compensation
            
            %N=20000;
            
            % First: find integer alignment
            r=xcorr(x,y);
            r2=sort(abs(r));
            if r2(floor(end*0.99))>0.3*r2(end),warning('The imput signals are not correlated, timealign cannot be performed'); end
            N=(length(r)-1)/2;
            [~,maxind]=max(abs(r));
            adjust=N+1-maxind;
            
            % Then: find accurate alignment to sub-sample resolution
            tdelta=fminbnd(@(tdelta) -abs(circdelay_local(r,tdelta-adjust,N+1)), -0.5, 0.5);  %,optimset('TolX',1e-12) if better accuracy needed
            
            % Do the time alignment
            y2=circdelay_local(y,adjust-tdelta);
            
            % Do the phase alignment
            L=min(length(x),length(y2));
            y2=y2.*exp(1i*angle(y2(1:L)'*x(1:L)));
            
            if isreal(y), y2=real(y2); end
            
            if nargout>1, timeadjustment=adjust-tdelta; end
            
            if nargin>2 && strcmpi(mode,'phasetracking')
                p=phasetracker(y2,x,1,1e-7);
                y2=y2.*exp(-1i*p);
            end
            
            function x2=circdelay_local(x,delay,N)
                % time shifting by a linear phase addition in the spectral domain.
                % Parameter N is only used in the optimization, otherwise unnessesary
                
                x2=ifft(ifftshift(fftshift(fft(x)).*exp(1i*2*pi*delay*(-length(x)/2:length(x)/2-1)'/length(x))));
                if nargin>2, x2=x2(N); end % this is just for the optimization fminbnd
                
            end
            
        end
        function [y, Pcon, Pout] = weblab(x, vg)
            
            if nargin<2
                vg=-2.3;
            end
            
            %[y,~,Idc,Udc]=RFWebLab_PA_meas_xs_VG(x, vg);
            RMSin=Usefulfunctions.getdbm(x);
            [y,~,Idc,Udc]=RFWebLab_PA_meas_Vg_v1_1(x, RMSin,vg);
            
            %            y=Usefulfunctions.timealign(x,y,'phasetracking');
            y=Usefulfunctions.timealign(x,y);
            y=double(y);
            
            Pcon=Udc*Idc;
            Pout=mean(abs(y).^2)/50/2;
            %gain=mean(abs(y))/mean(abs(u));
            %eff=Pout/Pcon;
        end
    end
    methods (Static) % Communication functions
        function x=encode_concurrent(X,fc,bw)
            % usage:
            % > X=randn_c(1000,3);X=resample(X,3,1);
            % > x=encode_concurrent(X,[-0.3 0.1 0.2],0.1);
            % > y=pa(x);
            % > X2=decode_concurrent(y,[-0.3 0.1 0.2],0.1);
            % > nmse(X,X2)
            
            assert(size(X,2)==length(fc),'Size mismatch');
            
            if bw<1
                X=sincresample([zeros(20,size(X,2));X;zeros(20,size(X,2))],1/bw);
            end
            X3=X.*exp(1i*2*pi*fc(:)'.*time(X));
            x=sum(X3,2);
        end
        function X=decode_concurrent(x,fc,bw)
            
            X=x.*exp(-1i*2*pi*fc(:)'.*time(x));
            
            if bw<1
                a=fir1(256,bw);
                X=filtfilt(a,1,X);
                X=sincresample(X,bw);
            end
            X=X(21:end-20,:); % remove zeros added in encode_concurrent
            
        end
        function out=evm(x,y)
            % evm(x,y)): Error Vector Magnitude
            % The output is given as a ratio value (low is good, 0 is perfect).
            % DC is removed and gain is optimized before calculation
            % Approximately 20 dB is needed for 64-QAM, then 3 dB up/down per
            % halving/doubling of the constellation.
            assert(size(x,2)==1,'wrong size');
            assert(size(y,2)==1,'wrong size');
            assert(size(x,1)==size(y,1),'wrong size');
            
            H=[ones(size(y)) y];
            theta=pinv(H)*x;
            %out=snr(x,H*theta);
            
            out=std(x-H*theta)/std(x);
        end
        function par=hardwareanalyzer(x,y,tocompensate)
            % par=hardwareanalyzer(x,y): Tool to analyze hardware problems (Gain
            % offset, DC offset, AM/AM and AM/PM nonlinearity and I/Q imbalance,
            % x is the transmitted IQ data vector, y is the received IQ data vector.
            % x and y must be of equal length and time synchronized.
            x=x(:);y=y(:);
            if nargin<3
                tocompensate='pid'; % compensate gain, phase, iq, dc offset
            end
            par=struct;
            par.tocompensate=tocompensate;
            
            C=y;
            if strfind(par.tocompensate,'i'), C=[C conj(y)]; end
            if strfind(par.tocompensate,'d'), C=[C ones(size(y))]; end
            if strfind(par.tocompensate,'3'), C=[C y.*abs(y).^2]; end
            if strfind(par.tocompensate,'5'), C=[C y.*abs(y).^4]; end
            if strfind(par.tocompensate,'f'), C=[C 1i*linspace(-1,1,length(y))'.*y (linspace(-1,1,length(y))').^2.*y]; end
            if strfind(par.tocompensate,'e') 
                for del=-16:16
                    C=[C Usefulfunctions.delay(y,del)];
                end
            end
            par.coeff=pinv(C)*x;
            
            par.importance=sqrt(mean(abs(C).^2).'.*abs(par.coeff).^2);
        end
        function y_compensated=hardwarecompensator(y,par)
            % y_compensated=hardwarecompensator(y,par): Tool to compensate hardware problems (Gain
            % offset, DC offset, AM/AM and AM/PM nonlinearity and I/Q imbalance, and frequency).
            % y is the received IQ data vector, par is compensation parameters, as previously
            % derived by the hardwareanalyzer function.
            
            C=y;
            
            if strfind(par.tocompensate,'i'), C=[C conj(y)]; end
            if strfind(par.tocompensate,'d'), C=[C ones(size(y))]; end
            if strfind(par.tocompensate,'3'), C=[C y.*abs(y).^2]; end
            if strfind(par.tocompensate,'5'), C=[C y.*abs(y).^4]; end
            if strfind(par.tocompensate,'f'), C=[C 1i*linspace(-1,1,length(y))'.*y (linspace(-1,1,length(y))').^2.*y]; end
            if strfind(par.tocompensate,'e')
                for del=-16:16
                    C=[C Usefulfunctions.delay(y,del)];
                end
            end
            y_compensated=C*par.coeff;
        end
        function H=Hscattering(pos_array, pos_scatterer, pos_ue)
            % See documentation in Hlos
            % all positions given in lambdas
            
            H = Hlos(pos_scatterer,pos_ue) * Hlos(pos_array,pos_scatterer);
        end
        function H=Hlos(pos_array, pos_ue)
            % H=Hlos(pos_array,pos_ue);
            % Computes the channel matrix H between an antenna array and a set of terminals
            % (ue). The positions of antennas and users are given as complex numbers,
            % representing the coordinates in a plane. (in number of lambdas).
            % Antenna gains G1 and G2 are supposed to be 0 (isotropic).
            %
            % Example use:
            % p1=1i*[0;0.5;1;1.5;2];
            % p2=[100;100+1i*10];
            % H=Hlos(p1,p2);
            % See also Hlos3D.
            assert(min(size(pos_array))==1,'wrong vector dimensions');
            pos_array=pos_array(:);
            pos_ue=pos_ue(:);
            dist=abs(repmat(pos_array.',length(pos_ue),1 )-  repmat(pos_ue,1,length(pos_array)) );
            H=(1/4/pi./dist).*exp(1i*2*pi*dist);
        end
        function H=Hlos3D(pos_array, pos_ue, Hplane, Eplane)
            % H=Hlos3D(pos_array,pos_ue, Hplane, Eplane);
            % Computes the channel matrix H between an antenna array and a set of terminals
            % (ue). The pos_array and pos_ue parameters should be Nx3-dimensional matrices, where the two first
            % columns is the placement in the antenna plane, and the third columm shows
            % the distance between the RX and TX antennas.
            % The Hplane and Eplane arguments can be omitted. They indicate
            % the antenna diagram, in dB, at different angles, in both the
            % horizontal and vertical plane. They should both be 2-D matrices,
            % with the 1st column indicating angles (in radians), and the 2nd column the
            % loss in dB at that angle.
            
            assert(size(pos_array,2)==3 && size(pos_ue,2)==3,'Wrong antenna position matrices');
            if nargin>2, assert(size(Hplane,2)==2 && size(Eplane,2)==2,'Wrong Hplane/Eplane matrices'); end
            A1=reshape(repmat(pos_array,1,size(pos_ue,1))',3,size(pos_array,1)*size(pos_ue,1))';
            B1=reshape(repmat(pos_ue,size(pos_array,1),1),size(pos_array,1)*size(pos_ue,1),3);
            dist=sqrt(reshape(sum((A1-B1).^2,2),size(pos_ue,1),size(pos_array,1)));
            H=(1/4/pi./dist).*exp(1i*2*pi*dist);
            
            % what angles do you need?
            fi_h=zeros(size(pos_ue,1),size(pos_array,1));
            fi_v=zeros(size(pos_ue,1),size(pos_array,1));
            for t=1:size(pos_array,1)
                for r=1:size(pos_ue,1)
                    dx=pos_ue(r,1)-pos_array(t,1);
                    dy=pos_ue(r,2)-pos_array(t,2);
                    dz=pos_ue(r,3)-pos_array(t,3);
                    fi_h(r,t)=atan(dx/dz);
                    fi_v(r,t)=atan(dy/dz);
                end
            end
            if nargin>2
                LossH=interp1(Hplane(:,1),Hplane(:,2),fi_h,'Cubic');
                LossV=interp1(Eplane(:,1),Eplane(:,2),fi_v,'Cubic');
                H=H.*10.^(LossH/20+LossV/20);
            end
        end
        function linkdoctor(x,y)
            % theta=linkdoctor(x,y)
            % The inputs x and y are baseband signals.
            % Normally they should be sampled at symbol rate, but the code works also
            % with oversampled baseband.
            % (c) Thomas Eriksson 2018
            % Modified:
            %   added bandwidth estimation (jan 2020)
            
            n=0;
            W=bandwidth_estimate_local(x); % BW estimate. Returns 1 if no oversampling.
            H=[];
            
            n=n+1;T{n}='Scaling                  ';
            H(:,n)=x;
            
            n=n+1;T{n}='DC offset                ';
            H(:,n)=ones(size(x));
            
            n=n+1;T{n}='I/Q imbalance           ';
            H(:,n)=conj(x);
            
            n=n+1;T{n}='Carrier frequency offset ';
            H(:,n)=1i*x.*linspace(-1,1,length(x))';
            
            n=n+1;T{n}='CFO and IQI              ';
            H(:,n)=1i*conj(x).*linspace(-1,1,length(x))';
            
            n=n+1;T{n}='Nonlinearity (3rd order) '; % should be matched filtered. Perhaps ideal LP?
            H(:,n)=x.*abs(x).^2;
            
            n=n+1;T{n}='Delay                    ';
            O=256;
            del=(-O:O)';
            h=(-1).^del./del;h(O+1)=0;
            h=h.*hanning(length(h));
            temp=filter(h,1,[x;zeros(O,1)]);
            temp=temp(O+1:end);
            H(:,n)=temp;
            
            n=n+1;T{n}='Differing delay I and Q  ';
            H(:,n)=conj(temp);
            
            % frequency-selective IQI
            
            n1=n;
            T{n+1}='Frequency selectivity    ';
            for i=-10:10  % frequency selectivity
                if i~=0
                    n=n+1;
                    if i<=0, H(:,n)=[x(-i+1:end);zeros(-i,1)];
                    else, H(:,n)=[zeros(i,1);x(1:end-i)]; end
                end
            end
            
            % phase noise missing
            % wrong with nonlinearities
            % complex conjugate CFO?
            %--------------------------------
            [Q,~]=qr(H,0);
            
            theta=Q'*y;
            
            %th=[ th(1:n1) ; sum(th(n1+1:end)) ]; % The effect of the frequency selectivity is summed.
            th=[ theta(1:n1) ; sqrt(abs(sum(theta(n1+1:end)).^2)) ]; % The effect of the frequency selectivity is summed.
            
            %[~,ind]=sort(-abs(th));
            
            disp('The following is a list of the dominating effects, given as EVM [%]');
            disp(['The relative bandwidth is ',num2str(W)]);
            disp('-------------------------------------------------------------------');
            
            N=abs(th)/abs(th(1));
            for i=2:length(th)
                disp([T{i},num2str(round(N(i)*100,1)),'%']);
            end
            
            e=y-Q*theta;
            H2=[ones(size(x)) abs(x) abs(x).^2 abs(x).^3];
            H2=H2./sqrt(sum(abs(H2).^2));
            [Q2,~]=qr(H2,0);
            theta2=pinv(Q2)*abs(e);
            
            D=sqrt(sum(abs(y-Q*theta).^2)/sum(abs(y).^2));
            D1=D*(abs(theta2(1)).^2/sum(abs(theta2).^2));
            D2=D*(abs(theta2(2)).^2/sum(abs(theta2).^2));
            D3=D*(abs(theta2(3)).^2/sum(abs(theta2).^2));
            D4=D*(abs(theta2(4)).^2/sum(abs(theta2).^2));
            
            disp(' ');
            disp(['Multiplicative noise |x|    ',num2str(round(D2*100,1)),'%']);
            disp(['Multiplicative noise |x|^2  ',num2str(round(D3*100,1)),'%']);
            disp(['Multiplicative noise |x|^3  ',num2str(round(D4*100,1)),'%']);
            disp(['Additive noise              ',num2str(round(D1*100,1)),'%']);
            function [W,fc,linf,SNR]=bandwidth_estimate_local(x)
                % Bandwidth wstimation. Returns 1 if the signal is using the full BW, <1
                % otherwise. The code can also estimate center frequency and approximate
                % noise level and SNR (if oversampled, such that parts of the spectrum is noise only).
                
                X=spec_local(x);
                Xl=10.^(spec_local(x)/10);
                f=(1:length(X))';f=f-mean(f);
                p=zeros(4,1);
                p(2)=f'*Xl/sum(Xl);
                m1=sum(Xl);
                m2=sum(Xl.^2);
                p(3)=m1^2/m2;
                
                
                H=[rect_local((f-p(2))/p(3)) ones(size(f))];
                th=pinv(H)*X;
                p(1)=th(1);
                p(4)=th(2);
                
                W=fminsearch(@(W) sum(abs(Xl-10.^( (p(1).*rect_local((f-p(2))./W)+p(4))/10  )).^2),p(3));
                popt=fminsearch(@(p) sum(abs(X-p(1).*rect_local((f-p(2))./W)-p(4)).^2),p);
                
                W=W/length(X);
                fc=popt(2)/length(X);
                linf=popt(4);
                SNR=popt(1);
                
                if SNR<3
                    W=1;
                    SNR=inf;
                    linf=-inf;
                end
            end
            function out=spec_local(inputsignal,N,fs)
                
                if nargin<3 fs=1; end
                if nargin<2 N=1024; end
                
                if size(inputsignal,2)>size(inputsignal,1)
                    inputsignal=inputsignal';
                end
                if (length(inputsignal)<N)
                    inputsignal=[zeros(floor((N-length(inputsignal))/2),size(inputsignal,2));inputsignal;zeros(ceil((N-length(inputsignal))/2),size(inputsignal,2))];
                end
                spectrumplotcolor=[0 0 1 ; 1 0 0 ; 0 1 0 ; 0 0 0 ; 1 1 0 ; 1 0 1; 0 1 1; 0.6 0.6 0.6];
                spectrumplotcolor=[spectrumplotcolor;0.6*spectrumplotcolor];
                legendstrings=[];
                for k=1:size(inputsignal,2)
                    x=inputsignal(:,k);
                    N=min(N,length(x)-1);
                    [X,w]=pwelch(x*sqrt(2*pi),hanning(N),[],'centered');
                    if nargout==1
                        out(:,k)=10*log10(X);
                    else
                        plot(w/2/pi*fs,10*log10(X),'color',spectrumplotcolor(k,:),'linewidth',2);
                        xlabel('Normalized frequency');
                        ylabel('Spectrum [dB]');
                        grid on
                        hold on
                    end
                    legendstrings=[legendstrings;num2str(k)];
                end
                if nargout~=1
                    axis tight;
                    legend(legendstrings);
                    hold off
                end
            end
            function out=rect_local(t)
                
                out=zeros(size(t));
                out(abs(t)<0.5)=1;
                
            end
            
        end
        function const=makespiralconstellation(M,f)
            % f=bestf(M,snr(i),var_pn)/2;
            
            m=(1:M)';
            ang=sqrt( (4*pi*m).^2*f/2+sqrt( ((4*pi*m).^2*f/2).^2+(4*pi*m).^2 ) );
            
            const=ang.*exp(1i*ang);
            
        end
        function xmf = matchedfilter(x,oversampling,beta,FL)
            % x is a vecor or matrix where the columns are samples to filter.
            % The resulting vector contains extra samples for filling and emtying the
            % filters in the pulseshape and matchedfilter functions. Remove
            % 2*oversampling*FL samples from the beginning of the vector
            % The matchedfilter function do not provide subsampling at the end
            
            if nargin<4, FL=25; end
            if nargin<3, beta=0.07; end
            if nargin<2, oversampling=5;end
            if FL<round(1/beta)
                disp('Warning: FilterLength too small to guarantee 35 dB');
            end
            
            PulseFilter=Usefulfunctions.rrc(0:1/oversampling:FL,beta);
            PulseFilter=[fliplr(PulseFilter(2:end)) PulseFilter]; % this makes sure that t=0 is represented
            
            xz = [ x ; zeros((length(PulseFilter)-1)/2,size(x,2)) ] ;
            
            xmf = filter(PulseFilter,1,xz)/oversampling;
        end
        function xps=pulseshape(x,oversampling,beta,FL,h)
            % pulseshape(x,oversampling,beta)
            % Root-raised-cosine pulse shaping.
            % x is a vecor or matrix where the columns are samples to pulseshape.
            % beta is the rootraisedcosine factor.
            % FL is the half filterlength, in symbols. Default 25.
            % h is pulse shape. Default rrc.
            %
            % The function oversamples and filters with rrc(beta) filter.
            % The resulting vector contains extra samples for filling and emtying the
            % filters in the pulseshape and matchedfilter functions. Remove
            % 2*oversampling*FL samples from the beginning of the vector
            
            if nargin<4, FL=25; end
            if nargin<3, beta=0.07; end
            if nargin<2, oversampling=5;end
            if FL<round(1/beta)
                disp('Warning: FilterLength too small to guarantee 35 dB.');
            end
            
            if Usefulfunctions.isint(oversampling)
                PulseFilter=Usefulfunctions.rrc(0:1/oversampling:FL,beta);
                PulseFilter=[fliplr(PulseFilter(2:end)) PulseFilter]; % this makes sure that t=0 is represented
                
                xz = [ x ; zeros(FL,size(x,2)) ] ; % append zeros for filter delay
                
                % multiphase filtering. Same as zeropadding and filtering.
                xps=zeros(size(xz,1)*oversampling,size(xz,2));
                for d=1:oversampling
                    xps(d:oversampling:end,:)=filter(PulseFilter(d:oversampling:end),1,xz);
                end
            else
                x=[zeros(FL,1);x];
                %    t=linspace(1,length(x)+1-1/oversampling-0.000001,round(length(x)*oversampling));
                t=1:1/oversampling:length(x)+1-0.001;
                x=[zeros(FL,1);x;zeros(FL+1,1)];
                t=t+FL;
                xps=zeros(length(t),size(x,2));
                p1=round(t-FL);
                p2=round(t+FL);
                
                for k=1:length(t)
                    t2=(p1(k):p2(k))-t(k)+eps;
                    xps(k,:)=( sin(pi*t2*(1-beta)) + 4*beta*t2.*cos(pi*t2*(1+beta)) )./(pi*t2.*(1-(4*beta*t2).^2))*x(p1(k):p2(k),:);
                end
                
            end
        end
        function [x,c] = randconst( rows, cols, M, type)
            % RANDCONST Complex constellation
            if nargin<4, type='QAM'; end
            if nargin<3, M=16; end
            if nargin<2
                if length(rows)==2
                    cols=rows(2);
                    rows=rows(1);
                end
            end
            
            switch type
                case 'QAM'
                    switch M
                        case 2
                            c=[-1;1];
                        case 8
                            c=[-3 -1 1 3;-3 -1 1 3]+1i*[ones(1,4); -ones(1,4)];
                            c=c(:);
                        case 32
                            [X,Y]=meshgrid(-5:2:5,-3:2:3);
                            c=X(:)+1i*Y(:);
                            c=[ (-3:2:3)'-1i*5 ; c ;  (-3:2:3)'+1i*5];
                        case 128
                            [X,Y]=meshgrid(-11:2:11,-7:2:7);
                            c=X(:)+1i*Y(:);
                            c=[ (-7:2:7)'-1i*11 ; (-7:2:7)'-1i*9 ; c ;  (-7:2:7)'+1i*9 ; (-7:2:7)'+1i*11 ];
                        case 512
                            [X,Y]=meshgrid(-15:2:15,-23:2:23);
                            [X2,Y2]=meshgrid(-3:2:3,-15:2:15);
                            c=[X(:)+1i*Y(:);X2(:)+1i*Y2(:)-20;X2(:)+1i*Y2(:)+20];
                            %                c=[ (-15:2:15)'-1i*11 ; (-15:2:15)'-1i*9 ; c ;  (-7:2:7)'+1i*9 ; (-7:2:7)'+1i*11 ];
                        otherwise
                            q=log2(M);
                            if round(q/2)~=q/2
                                error('bad constellation size');
                            end
                            q=round(sqrt(M));
                            r=((1:q)-(q+1)/2);
                            c=reshape(repmat(r,q,1)+1i*repmat(r',1,q),q^2,1);
                    end
                case 'PSK'
                    c=exp(1i*2*pi*(1:M)/M).';
                case 'SPIRAL'
                    c=makespiralconstellation(M,0);
                otherwise
                    error('Wrong constellation type. Choose QAM, PSK or SPIRAL');
            end
            c=sqrt(2)/std(c,1)*c;
            
            
            i=randi(size(c,1),rows,cols);
            x=c(i);
        end
        function pulse=rc(t,beta)
            
            %i2=find(1-4*beta^2*t.^2<1e-13);t(i2)=0.11;
            i2=find((4*beta^2*t).^2==1);t(i2)=0.11;
            
            pulse=sinc(t).*cos(pi*beta*t)./(1-4*beta^2*t.^2);
            
            pulse(i2)=0;
        end
        function pulse=rrc(t,beta)
            
            i1=find(t==0);t(i1)=0.11; % value 0.11 not used
            i2=find((4*beta*t).^2==1);t(i2)=0.11;
            pulse=( sin(pi*t*(1-beta)) + 4*beta*t.*cos(pi*t*(1+beta)) )./(pi*t.*(1-(4*beta*t).^2));
            pulse(i1)=(1-beta+4*beta/pi);
            pulse(i2)=beta/sqrt(2)*((1+2/pi)*sin(pi/4/beta)+(1-2/pi)*cos(pi/4/beta));
            
            
        end
    end
    methods (Static) % Oscillator functions
        function xs=phasetracker(y, x, var_awgn, var_pn,a)
            % phase tracker program
            % Model: y = x * exp(i * xs) + w
            % y is the phase-noisy signal, to be cleaned
            % x is the (estimated) input signal
            % var_awgn is the variance of the error of y. This should be a vector of
            % the same length as y. For pilots, the var_awgn variable should be equal
            % to the channel variance, for data it should indicate how much the data
            % can be trusted (a value maybe 2-3 times the channel variance, but that
            % varies)
            % var_pn is the variance of the phase noise delta_fi.
            % The parameter a can be omitted (default value 1). It is the forgetting
            % factor in the pn filter. For Wiener, a=1.
            
            if nargin<5
                a=1;
            end
            y=y./x;                 % remove estimated symbol from y
            R=var_awgn/2./abs(x).^2;  % compensate variance accordingly
            
            % N=size(y,2);
            N=length(y);
            
            Pf=zeros(size(y));Pp=Pf;Ps=Pf;  % Initialize variables for Kalman filtering
            xf=zeros(size(y));xp=xf;xs=xf;
            
            
            % Kalman filtering
            xp(1)=angle(y(1));
            Pp(1)=real(1+var_pn);
            for i=1:N
                if i>1
                    xp(i)=a*xf(i-1);
                    Pp(i)=real( a^2*Pf(i-1)+var_pn);
                end
                H=1i*exp(1i*xp(i));
                K=Pp(i)*conj(H)/(Pp(i)+R(i));
                xf(i)=xp(i)+real(K*(y(i) - exp(1i*xp(i))));
                Pf(i)=real(Pp(i)-K*H*Pp(i));
            end
            
            % Kalman smoother; run it all backwards.
            xs(end)=xf(end);
            Ps(end)=Pf(end);
            
            for i=N-1:-1:1
                C=real(a*Pf(i)/Pp(i+1));
                xs(i)=real(xf(i)+C*(xs(i+1)-xp(i+1)));
                %    Ps(i)=real(Pf(i)+C^2*(Ps(i+1)-Pp(i+1)));
            end
            
        end
        function xs=phasortracker(y, x, var_awgn, var_pn)
            % phase tracker program
            % Model: y = x * exp(i * xs) + w
            % y is the phase-noisy signal, to be cleaned
            % x is the (estimated) input signal
            % var_awgn is the variance of the error of y. This should be a vector of
            % the same length as y. For pilots, the var_awgn variable should be equal
            % to the channel variance, for data it should indicate how much the data
            % can be trusted (a value maybe 2-3 times the channel variance, but that
            % varies)
            % var_pn is the variance of the phase noise delta_fi
            
            y=y./x;                 % remove estimated symbol from y
            R=var_awgn/2./abs(x).^2;  % compensate variance accordingly
            
            % N=size(y,2);
            N=length(y);
            
            Pf=zeros(size(y));Pp=Pf;Ps=Pf;  % Initialize variables for Kalman filtering
            xf=zeros(size(y));xp=xf;xs=xf;
            
            % Kalman filtering
            xp(1)=y(1);
            Pp(1)=real(1+var_pn);
            for i=1:N
                if i>1
                    xp(i)=xf(i-1);
                    Pp(i)=real(Pf(i-1)+var_pn);
                end
                K=Pp(i)/(Pp(i)+R(i));
                xf(i)=xp(i)+(K*(y(i) - xp(i)));
                Pf(i)=real(Pp(i)-K*Pp(i));
            end
            
            % Kalman smoother; run it all backwards.
            xs(end)=xf(end);
            Ps(end)=Pf(end);
            for i=N-1:-1:1
                xs(i)=xf(i)+real(Pf(i)/Pp(i+1))*(xs(i+1)-xp(i+1));
            end
            xs=angle(xs);
            
        end
        function [fout,Xout]=pnspec(x,fs)
            % (c) Thomas Eriksson 2018
            
            if ~isreal(x)
                x=x/sqrt(mean(abs(x).^2));
            end
            
            % compute power at 100 kHz, by filtering
            BW=20e3;
            x2=Usefulfunctions.fftbpfilter(x,100e3/fs,BW/fs);
            %Pdb=pow2db(var(x2)/0.3754/BW/2);
            Pdb=pow2db(var(x2)/BW/2);
            
            %'check this'
            disp(['Power at 100 kHz: ',num2str(Pdb)]);
            
            
            
            d=round((1:length(x)/50:length(x)/2)-1); % 50 steps
            f=[];SX=[];
            
            %SX2=(abs(fft(x)).^2)/length(x)/2/fs; % worst row
            
            
            % technique to create a close grid of values
            % by computing different-size fft's
            for i=1:length(d)
                x2=x(1:end-d(i));
                x3=x(1+d(i):end);
                SX2=(abs(fft(x2)).^2+abs(fft(x3)).^2)/length(x2)/2/fs; % worst row
                %    SX2=abs(fft(x)).^2/length(x)/2/fs; % worst row
                SX2=SX2(1:floor(length(SX2)/2));
                
                f2=linspace(0,0.5*fs,length(SX2));
                f=[f;f2'];
                SX=[SX;SX2];
            end
            
            % SX=abs(fft(x)).^2/length(x)/2/fs; % worst row
            % SX=SX(1:floor(length(SX)/2));
            % f=linspace(0,0.5*fs,length(SX));
            
            [~,k]=sort(f);
            f=f(k);
            SX=SX(k);
            
            
            
            % remove 0 frequency part
            firstind=find(f>0);
            f=f(firstind(1):end);
            SX=SX(firstind(1):end);
            
            
            
            steps=100;
            stepsize=logspace(-5,0,steps)';
            stepsize=max(1,round((stepsize/sum(stepsize)*length(SX))));
            stepsize(end)=stepsize(end)+(length(SX)-sum(stepsize));
            assert(sum(stepsize)==length(SX));
            
            SX2=zeros(steps,1);
            f2=zeros(steps,1);
            n=0;
            vv=[];
            vf=[];
            for i=1:steps
                S=SX(n+1:min(length(SX),n+stepsize(i)));
                SX2(i)=mean(S); % gives power per hertz. makes sense if spec is almost flat in region.
                f2(i)=mean(f(n+1:min(length(f),n+stepsize(i))));
                
                Ss=sort(S);
                PAPR=10*log10(abs(Ss(end))/abs(Ss(round(end*0.8))));
                if PAPR>10 && pow2db(mean(S))>-150 && pow2db(SX2(i)/SX2(i-1))>2
                    vv=[vv sum(S)*fs/length(SX)];
                    [~,ind]=max(S);
                    vf=[vf f(n+ind)];
                end
                n=n+stepsize(i);
            end
            
            semilogx(f2,pow2db(SX2),'-',1e5,Pdb,'k*',vf,pow2db(vv),'rx');
            grid
            
            % Add lines for spurioses
            mindb=min(pow2db(SX2));
            for i=1:length(vv)
                line([vf(i) vf(i)],[mindb pow2db(vv(i))],'color','r');
            end
            
            % Add a line indicating -20 dB/decade
            line([1e5/10^0.25 1e5*10^0.25],[Pdb+5,Pdb-5],'color','g')
            
            if nargout>0
                fout=f2;
                Xout=SX2;
            end
            
            M1=min(pow2db(SX2));
            M2=max(pow2db(SX2));
            Psum=cumsum(stepsize.*SX2);
            Psum=Psum/Psum(end)*(M2-M1)+M1;
            hold on
            semilogx(f2,Psum,'color',[0.8 0.8 0.8]);
            hold off
            
            
        end
        function pn(varargin)
            
            persistent l0 l100 linf W
            persistent hmain f3db htext100 htextinf htext3db hspec hsnr
            persistent editW checkkalman
            persistent l0button
            persistent x
            
            if nargin==0
                pn('initialize');
            else
                switch varargin{1}
                    case 'initialize'
                        clf;
                        linf=-140;
                        l100=-90;
                        f3db=10;
                        W=2e7;
                        figure(1);
                        hspec=axes('OuterPosition',[0 0.1 0.5 0.9]);
                        hsnr=axes('OuterPosition',[0.5 0.1 0.5 0.9]);
                        l0button = uicontrol('Style', 'pushbutton','String', 'quit',...
                            'units','normalized','position',[0 0.0 0.05 0.05],...
                            'backgroundcolor',[0.9 0.5 0.5],...
                            'Callback', 'pn( ''quit'' );' ,'visible','on');
                        
                        uicontrol('style','text','string','Sampling frequency [MHz] ','position',[80,20,130,20],'HorizontalAlignment','right');
                        checkkalman=uicontrol('style','checkbox','string','Kalman ','position',[500,10,90,20],'HorizontalAlignment','right','callback','pn(''edited'',1)');
                        
                        
                        
                        editW=uicontrol('style','edit','string',num2str(floor(W/1e5)/10),...
                            'position',[210,22,30,20],'callback','pn(''edited'',1)');
                        
                        pn('draw');
                        pn('drawsnr');
                        return;
                    case 'draw'
                        f=logspace(0,8);
                        
                        % uträkning av spectrum
                        l0=100+l100-20*log10(f3db);
                        S=db2pow(l0)*f3db^2./(f3db^2+f.^2)+db2pow(linf);
                        
                        % rita spektrum
                        axes(hspec);
                        hmain=semilogx(f,10*log10(S),'linewidth',5);
                        axis([f(1),f(end),-150,0]);
                        line([f(1),f(end)],[linf,linf],'color','r','linestyle',':');
                        line([f3db,f3db],[-150,0],'color','k','linestyle',':');
                        
                        % rita greppunkter
                        %line([1e5/1.25,1e5*1.25],[l100,l100],'linewidth',30,'ButtonDownFcn', 'pn(''mousepress'',1) ');
                        %htext100=text(1.5e5, l100+5, [num2str(l100) ' dBc@100 kHz']);
                        %line([1e3/1.25,1e3*1.25],[linf,linf],'linewidth',10,'color','r','ButtonDownFcn', 'pn(''mousepress'',2) ');
                        %htextinf=text(1.5e3, linf+5, [num2str(linf) '']);
                        %line([f3db/1.25,f3db*1.25],[-60,-60],'linewidth',10,'color','k','ButtonDownFcn', 'pn(''mousepress'',3) ');
                        %htext3db=text(1.5*f3db, -60+5, [num2str(f3db) '']);
                        
                        patch([1e5/1.4,1e5/1.4,1e5*1.4,1e5*1.4],[l100-3,l100+3,l100+3,l100-3],[0 0 0],'ButtonDownFcn', 'pn(''mousepress'',1) ');
                        htext100=text(1.5e5, l100+5, [num2str(l100) ' dBc@100 kHz']);
                        
                        patch([1e3/1.4,1e3/1.4,1e3*1.4,1e3*1.4],[linf-3,linf+3,linf+3,linf-3],[1 0 0],'ButtonDownFcn', 'pn(''mousepress'',2) ');
                        htextinf=text(1.5e3, linf+5, [num2str(linf) ' dBc']);
                        
                        patch([f3db/1.4,f3db/1.4,f3db*1.4,f3db*1.4],[-60-3,-60+3,-60+3,-60-3],[0 0 0],'ButtonDownFcn', 'pn(''mousepress'',3) ');
                        htext3db=text(1.5*f3db, -60+5, [num2str(f3db) ' Hz']);
                        
                        
                        return;
                    case 'drawsnr'
                        
                        sigmadelta2=(1-exp(-4*pi*f3db/W))*pi*1e10*db2pow(l100)/f3db;
                        10*log10(sigmadelta2);
                        a=exp(-2*pi*f3db/W);
                        sn=10:3:40;
                        sigmaw2=10.^(-sn/10);
                        sigmae2=pi*1e5/sqrt(2)*sqrt(db2pow(l100))/sqrt(W).*sqrt(sigmaw2);
                        sigmae2_2=2*1e10*db2pow(l100)*sigmaw2.^2/W^2./(sigmaw2/W+2*db2pow(linf)).^1.5./sqrt(f3db^2*sigmaw2/W+2*db2pow(linf)*f3db^2+2*1e10*db2pow(l100)).*atan(W/2./sqrt(f3db^2+2e10*db2pow(l100)./(sigmaw2/W+2*db2pow(linf))))  +  sigmaw2.*db2pow(linf)./(sigmaw2/W+2*db2pow(linf));
                        axes(hsnr);
                        if get(checkkalman,'value')
                            evm_1=zeros(size(sn));
                            evm_2=zeros(size(sn));
                            for i=1:length(sn)
                                [a1,a2]=evaluate_pn(sn(i)-10*log10(W/20e6),sigmadelta2,20,4);
                                evm_1(i)=var(a1);
                                evm_2(i)=var(a2);
                            end
                            plot(sn,-10*log10(evm_1),'k:',sn,-10*log10(evm_2),'k');
                            hold on
                        end
                        plot(sn,-10*log10(sigmae2_2),'r:',sn,-10*log10(sigmae2),'r');
                        hold off;
                        axis([sn(1) sn(end) 20 60])
                        grid on
                        sn(end);
                        -10*log10(sigmae2_2(end));
                        return;
                        
                    case 'mousepress'
                        selectiontype=get(gcf,'selectiontype');
                        set(gcf,'windowbuttonmotionfcn',['pn(''mousemove'',' num2str(varargin{2}), ') ']);
                        set(gcf,'windowbuttonupfcn',['pn(''mouseup'',' num2str(varargin{2}), ') ']);
                        return
                    case 'mousemove'
                        set(gcf,'windowbuttonmotionfcn','');
                        pt=get(hspec,'currentpoint');
                        switch varargin{2}
                            case 1
                                l100=min(-40,max(linf,round(pt(1,2))));
                            case 2
                                linf=max(-150,min(l100,round(pt(1,2))));
                            case 3
                                f3db=round(pt(1,1));
                        end
                        pn('draw');
                        if ~get(checkkalman,'value')
                            pn('drawsnr');
                        end
                        set(gcf,'windowbuttonmotionfcn',['pn(''mousemove'',' num2str(varargin{2}), ') ']);
                        return
                    case 'mouseup'
                        set(gcf,'windowbuttonmotionfcn','');
                        set(gcf,'windowbuttonupfcn','');
                        pn('drawsnr');
                        return
                    case 'edited'
                        W=str2double(get(editW,'string'))*1e6;
                        pn('draw');
                        pn('drawsnr');
                    case 'quit'
                        close(1);
                        return;
                end
            end
            
            function [a1,a2]=evaluate_pn(SNR,var_pn,pilot_interval,M,a,METHOD)
                % EVALUATE evaluates the perfromance of phase tracking in a communication
                % system.
                % EVALUATE(SNR,VAR_PN,PILOT_INTERVAL) evaluates the performance with
                % given settings of the AWGN SNR, phase noise variance, and the interval
                % which to send pilots.
                % Example: Evaluate the performance for a phase tracker at 30 dB SNR, PN
                % variance 1e-4 rad^2, pilots every 20 samples.
                % e=evaluate(30, 1e-4, 20);
                % The output a1 and a2 are the entire estimation error sequences, before and after decision feedback.
                % The variance of the remaining phase error can thus be computed by
                % var(a2).
                
                if nargin<6
                    METHOD=1;
                    if nargin<5
                        a=1;
                        if nargin<4
                            M=16;
                            if nargin<3
                                pilot_interval=20;
                            end
                        end
                    end
                end
                DECISIONFEEDBACK=1;
                
                % Initializations
                S = RandStream.getGlobalStream;S.reset;
                
                h=modem.qammod('M',M);
                g=modem.qamdemod(h);
                
                var_signal=var(h.constellation,1);
                var_awgn=10^(-SNR/10)*var_signal;
                
                % pilots
                N=100000; % set approx nuber of iterations
                
                pilot_seq=1:pilot_interval:N;
                N=pilot_seq(end); % Adjust N so that the last symbol is a pilot
                
                % encoder
                index=randi([0,M-1],1,N);
                index(pilot_seq)=0; % ensures max power of pilots
                x=modulate(h,index);
                
                % channel
                pn_r=exp(1i*(cumsum(sqrt(var_pn)*randn(1,N),2) + repmat(2*pi*rand(1,1),1,N)));
                awgn=sqrt(var_awgn/2)*(randn(1,N)+1i*randn(1,N)); %AWGN
                
                y=x.*pn_r+awgn;
                
                % display true variance
                % disp(['True SNR: ' num2str((10*log10((sum(y.*conj(y),2)-sum(awgn.*conj(awgn),2))./sum(awgn.*conj(awgn),2)))')]);
                
                % phase follower
                var_awgn_sequence=1e99*ones(size(y));var_awgn_sequence(pilot_seq)=var_awgn; % Trust pilots up to awgn
                xhat=ones(size(y));xhat(pilot_seq)=h.constellation(1);
                switch METHOD,
                    case 0, phase_hat=phasetracker(y,xhat,var_awgn_sequence,var_pn,a);
                    case 1, phase_hat=phasortracker(y,xhat,var_awgn_sequence,var_pn);
                    case 2, phase_hat=MAP_PN_Estimator(y.',xhat.',var_awgn_sequence.',var_pn).';
                end
                xhat=y.*exp(-1i*phase_hat);
                
                indexhat=demodulate(g,xhat);
                indexhat(pilot_seq)=0;
                
                % t=0:length(C):prod(size(D))-1;
                % D(indexhat+t)=1e99;
                % [~,indexhat]=min(D);
                % D(indexhat+t)=1e99;
                % [~,indexhat]=min(D);
                % D(indexhat+t)=1e99;
                % [~,indexhat]=min(D);
                
                a1=(wrap(angle(pn_r)-phase_hat));
                
                % Hard decision feedback
                
                if DECISIONFEEDBACK==1
                    xhat2=modulate(h,indexhat);
                    var_awgn_sequence=3*var_awgn*ones(size(y));var_awgn_sequence(pilot_seq)=var_awgn;
                    switch METHOD,
                        case 0, phase_hat=phasetracker(y,xhat2,var_awgn_sequence,var_pn,a);
                        case 1, phase_hat=phasortracker(y,xhat2,var_awgn_sequence,var_pn);
                        case 2, phase_hat=MAP_PN_Estimator(y.',xhat2.',var_awgn_sequence.',var_pn).';
                    end
                    xhat2=y.*exp(-1i*phase_hat);
                    
                    indexhat2=demodulate(g,xhat2);indexhat2(pilot_seq)=0;
                    
                    EVM1=round(std(x-xhat)/std(x)*1000)/10;
                    EVM2=round(std(x-xhat2)/std(x)*1000)/10;
                    
                    SER1=sum(index~=indexhat)/length(index);
                    SER2=sum(index~=indexhat2)/length(index);
                    %    disp(['Error vector magnitude (EVM), only pilots used for tracking: ' num2str(EVM1) '%   EVM using decision feedback: ' num2str(EVM2) '%']);
                    %    disp(['Symbol error rate (SER): ' num2str(SER1) '   SER using decision feedback: ' num2str(SER2)]);
                end
                
                % performance
                
                %disp(['Raw symbol error rate : ' num2str(SER)]);
                %disp(['SNR after phase compensation: ' num2str(10*log10(mean(abs(x(pilot_seq_inv)).^2)/mean(abs(x(pilot_seq_inv)-xhat(pilot_seq_inv)).^2)))]);
                
                a2=(wrap(angle(pn_r)-phase_hat));
                
                
                
                function y  = wrap(x)
                    %#eml
                    f=2*pi;
                    t=floor((x+pi)/f);
                    y=x-f*t;
                    
                end
            end
        end
        function pn2(varargin)
            
            persistent freq level spurfreq spurlevel W
            persistent hmain hspec
            persistent editW buttonSpur editSpur
            %persistent l0button htext
            
            if nargin==0
                pn2('initialize');
            else
                switch varargin{1}
                    case 'initialize'
                        clf;
                        freq=[1e3 1e4 1e5 1e6];
                        level=[-40 -80 -100 -120];
                        spurfreq=[1e5 2e5 3e5];
                        spurlevel=[-20 -30 -40];
                        
                        W=2e7;
                        figure(1);
                        hspec=axes('OuterPosition',[0 0.1 0.9 0.9]);
                        l0button = uicontrol('Style', 'pushbutton','String', 'quit',...
                            'units','normalized','position',[0 0.0 0.05 0.05],...
                            'backgroundcolor',[0.9 0.5 0.5],...
                            'Callback', 'pn2( ''quit'' );' ,'visible','on');
                        
                        uicontrol('style','text','string','Max freq [MHz]   ','position',[80,20,130,20],'HorizontalAlignment','right');
                        
                        editW=uicontrol('style','edit','string',num2str(floor(W/1e5)/10),...
                            'position',[210,22,30,20],'callback','pn2(''edited'',1)');
                        
                        
                        buttonSpur=uicontrol('style','checkbox','string','Add spurs at',...
                            'position',[270,22,80,20],'callback','pn2(''pressed'',1)');
                        editSpur=uicontrol('style','edit','string',num2str(floor(spurfreq(1)/1e3)),...
                            'position',[350,22,40,20],'callback','pn2(''spuredited'',1)');
                        uicontrol('style','text','string','kHz','position',[393,18,40,20],'HorizontalAlignment','left');
                        
                        uicontrol('Style', 'pushbutton','String', 'Generate PN',...
                            'position',[430,22,80,20],'backgroundcolor',[0.5 0.9 0.5],'visible','on');
                        uicontrol('Style', 'pushbutton','String', 'Plot BER',...
                            'position',[520,22,80,20],'backgroundcolor',[0.5 0.5 0.9],'visible','on');
                        
                        pn2('draw');
                        return;
                    case 'draw'
                        f=logspace(0,8,100);
                        
                        % uträkning av spectrum
                        S=interp1(log([1,freq,f(end)]),[level(1),level,level(end)],log(f),'pchip');
                        
                        
                        
                        % rita spektrum
                        axes(hspec);
                        hmain=semilogx(f,S,'linewidth',2);
                        grid on
                        axis([f(1),f(end),-150,0]);
                        
                        
                        f2=[f spurfreq spurfreq+1];
                        S2=[S spurlevel spurlevel];
                        [~,ind]=sort(f2);
                        f2=f2(ind);
                        S2=S2(ind);
                        M1=min(S2);
                        M2=max(S2);
                        stepsize=diff(f2);stepsize=[stepsize(1) stepsize];
                        Psum=cumsum(stepsize.*db2pow(S2));
                        Psum=Psum/Psum(end)*(M2-M1)+M1;
                        hold on
                        semilogx(f2,Psum,'color',[0.6 0.6 0.6]);
                        hold off
                        
                        
                        
                        
                        for i=1:length(freq)
                            patch(freq(i)*[1/1.4,1/1.4,1.4,1.4],level(i)+[-3,+3,+3,-3],[0 0 0],'ButtonDownFcn', ['pn2(''change_freqlevel'',',num2str(i),') ']);
                        end
                        
                        Sspurs=interp1(log([1,freq,f(end)]),[level(1),level,level(end)],log(spurfreq),'pchip');
                        for i=1:length(spurfreq)
                            line([spurfreq(i) spurfreq(i)],[Sspurs(i) spurlevel(i)],'color',[1 0 0],'linewidth',2);
                            patch(spurfreq(i)*[1/1.4,1/1.4,1.4,1.4],spurlevel(i)+[-3,+3,+3,-3],[1 0 0],'ButtonDownFcn', ['pn2(''change_spurfreqlevel'',',num2str(i),') ']);
                        end
                        
                        
                        set(editSpur,'string',num2str(floor(spurfreq(1)/1e3)));
                        set(gca,'ButtonDownFcn','pn2(''add'')');
                        set(hmain,'ButtonDownFcn','pn2(''add'')');
                        return;
                        
                    case 'change_freqlevel'
                        selectiontype=get(gcf,'selectiontype');
                        index=varargin{2};
                        if strcmp(selectiontype,'alt')
                            if length(freq)>=1
                                freq=freq([1:index-1,index+1:end]);
                                level=level([1:index-1,index+1:end]);
                            end
                            pn2('draw');
                        else
                            set(gcf,'windowbuttonmotionfcn',['pn2(''mousemove'',' num2str(index), ') ']);
                            set(gcf,'windowbuttonupfcn',['pn2(''mouseup'',' num2str(index), ') ']);
                        end
                        return
                        
                        
                    case 'change_spurfreqlevel'
                        selectiontype=get(gcf,'selectiontype');
                        index=varargin{2};
                        if strcmp(selectiontype,'alt')
                            spurfreq=spurfreq([1:index-1,index+1:end]);
                            spurlevel=spurlevel([1:index-1,index+1:end]);
                            pn2('draw');
                        else
                            set(gcf,'windowbuttonmotionfcn',['pn2(''mousemove2'',' num2str(index), ') ']);
                            set(gcf,'windowbuttonupfcn',['pn2(''mouseup'',' num2str(index), ') ']);
                        end
                        return
                        
                    case 'mousemove'
                        set(gcf,'windowbuttonmotionfcn','');
                        pt=get(hspec,'currentpoint');
                        
                        newf=min(max(round(pt(1,1)),1.1),1e8/1.1);
                        newl=round(pt(1,2));
                        ind=varargin{2};
                        if ind<length(freq) newf=min(newf,freq(ind+1)/1.1); end
                        if ind>1 newf=max(newf,freq(ind-1)*1.1); end
                        freq(ind)=newf;
                        level(ind)=newl;
                        
                        
                        
                        pn2('draw');
                        set(gcf,'windowbuttonmotionfcn',['pn2(''mousemove'',' num2str(varargin{2}), ') ']);
                        return
                        
                        
                    case 'mousemove2'
                        set(gcf,'windowbuttonmotionfcn','');
                        pt=get(hspec,'currentpoint');
                        
                        newf=min(max(round(pt(1,1)),1.1),1e8/1.1);
                        newl=round(pt(1,2));
                        ind=varargin{2};
                        if ind==1
                            spurfreq(1)=newf;
                            spurfreq(2)=newf*2;
                            spurfreq(3)=newf*3;
                        end
                        spurlevel(ind)=newl;
                        
                        
                        
                        pn2('draw');
                        set(gcf,'windowbuttonmotionfcn',['pn2(''mousemove2'',' num2str(varargin{2}), ') ']);
                        return
                        
                    case 'mouseup'
                        set(gcf,'windowbuttonmotionfcn','');
                        set(gcf,'windowbuttonupfcn','');
                        return
                    case 'add'
                        selectiontype=get(gcf,'selectiontype');
                        if strcmp(selectiontype,'alt')
                            pt=get(hspec,'currentpoint');
                            newf=round(pt(1,1));
                            newl=round(pt(1,2));
                            
                            f=logspace(0,8);
                            newl2=interp1(log([1,freq,f(end)]),[level(1),level,level(end)],log(newf),'pchip');
                            if abs(newl-newl2)<10
                                
                                freq=[freq newf];
                                level=[level newl];
                                [~,ind]=sort(freq);
                                freq=freq(ind);
                                level=level(ind);
                                pn2('draw');
                            end
                        end
                        return
                    case 'edited'
                        W=str2double(get(editW,'string'))*1e6;
                        pn2('draw');
                    case 'spuredited'
                        newf=str2double(get(editSpur,'string'))*1e3;
                        spurfreq=[newf,2*newf,3*newf];
                        pn2('draw');
                    case 'quit'
                        close(1);
                        return;
                end
            end
            
        end
    end
end