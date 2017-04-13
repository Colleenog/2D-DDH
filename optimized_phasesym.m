
function [OPS] = optimized_phasesym(im,nscale,norient,angle,minW,sigmaOnf, dThetaOnSigma,k,polarity);


% im= double(imageUS);nscale=1;norient=3;minWaveLength=25;mult=1;sigmaOnf=.25;
% dThetaOnSigma=2;k=3.5; polarity=1;orientWrap=0;thetaPhase = 1; .

im=double(im);
    v = version; Octave = v(1)<'5';  % Crude Octave test    
    epsilon         = .0001;         % Used to prevent division by zero.
 
    % Calculate the standard deviation of the angular Gaussian function
    % used to construct filters in the frequency plane.     
     thetaSigma = pi/norient/dThetaOnSigma;  
 
    % Ensure the number of rows and columns of im are even to ensure stable
    % behaviour with fftshift (if dimensions are odd fftshift(fftshift(a)) ~=  a).
    % Here I have simply clipped off one row or column as needed.
    
    [rows,cols] = size(im);
    crop = 0;
    if mod(rows,2)
    rows = rows-1;
    crop = 1;
    end
    if mod(cols,2)
    cols = cols-1;
    crop = 1;
    end
    if crop
    im = im(1:rows,1:cols);
    end
 
    imagefft = fft2(im);                % Fourier transform of image
    zero = zeros(rows,cols);
    
    totalEnergy = zero;                 % Matrix for accumulating weighted phase 
                                        % congruency values (energy).
    totalSumAn  = zero;                 % Matrix for accumulating filter response
                                        % amplitude values.
    orientation = zero;                 % Matrix storing orientation with greatest
                                        % energy for each pixel.
    estMeanE2n = [];
    EO = cell(nscale, norient);         % Cell array of convolution results
    ifftFilterArray = cell(1, nscale);  % Cell array of inverse FFTs of filters
 
    
    % Pre-compute some stuff to speed up filter construction
 
    [x,y] = meshgrid( [-cols/2:(cols/2-1)]/cols,...
                      [-rows/2:(rows/2-1)]/rows);
                  
    radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
    radius(rows/2+1, cols/2+1) = 1;   % Get rid of the 0 radius value in the middle 
                                      % so that taking the log of the radius will 
                                      % not cause trouble.
    theta = atan2(-y,x);              % Matrix values contain polar angle.
                                      % (note -ve y is used to give +ve
                                      % anti-clockwise angles)
    radius = fftshift(radius);        % Quadrant shift radius and theta so that filters
    theta  = fftshift(theta);         % are constructed with 0 frequency at the corners.
 
    sintheta = sin(theta);
    costheta = cos(theta);
    clear x; clear y; clear theta;    % save a little memory
 
    % Filters are constructed in terms of two components.
    % 1) The radial component, which controls the frequency band that the filter
    %    responds to
    % 2) The angular component, which controls the orientation that the filter
    %    responds to.
    % The two components are multiplied together to construct the overall filter.
    
    % Construct the radial filter components...
    
    % First construct a low-pass filter that is as large as possible, yet falls
    % away to zero at the boundaries.  All log Gabor filters are multiplied by
    % this to ensure no extra frequencies at the 'corners' of the FFT are
    % incorporated as this seems to upset the normalisation process when
    % calculating phase congrunecy.
    lp = lowpassfilter([rows,cols],.4,10);   % Radius .4, 'sharpness' 10
 
    logGabor = cell(1,nscale);
    
    logGabor = cell(nscale,norient);

    for p=1:norient
    for s = 1:nscale
        wavelength = minW(p)+5*(s-1);
        fo = 1.0/wavelength;                  % Centre frequency of filter.
        logGabor{s,p} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
        logGabor{s,p} = logGabor{s,p}.*lp;        % Apply low-pass filter
        logGabor{s,p}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                              % back to zero (undo the radius fudge).
    end
    end
 
    % Then construct the angular filter components...
    spread = cell(1,norient);
    

 
     angl=[pi/(180/angle(1)) pi/(180/angle(2)) pi/(180/angle(3))];
 
    for o = 1:norient
%         angl = (o-1)*pi/norient;           % Filter angle.
        
        % For each point in the filter matrix calculate the angular distance from
        % the specified filter orientation.  To overcome the angular wrap-around
        % problem sine difference and cosine difference values are first computed
        % and then the atan2 function is used to determine angular distance.
        
%         thetaSigma = pi/norient/dThetaOnSigma;
        ds = sintheta * cos(angl(o)) - costheta * sin(angl(o));    % Difference in sine.
        dc = costheta * cos(angl(o)) + sintheta * sin(angl(o));    % Difference in cosine.
        dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
        spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  % Calculate the
                                                             % angular filter component.
    end
 
    % The main loop...
 
    for o = 1:norient,                   % For each orientation.
        %fprintf('Processing orientation %d \r', o); 
        if Octave fflush(1); end
    
        sumAn_ThisOrient  = zero;      
        Energy_ThisOrient = zero;      
 
        for s = 1:nscale,                  % For each scale.
 
            filter = logGabor{s,o} .* spread{o};  % Multiply radial and angular
                                                % components to get filter.
 
            ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
            ifftFilterArray{s} = ifftFilt;                   % record ifft2 of filter
 
            % Convolve image with even and odd filters returning the result in EO
            EO{s,o} = ifft2(imagefft .* filter);
            An = abs(EO{s,o});                        % Amplitude of even & odd filter response.
            sumAn_ThisOrient = sumAn_ThisOrient + An; % Sum of amplitude responses.
            
            if s==1
                EM_n = sum(sum(filter.^2)); % Record mean squared filter value at smallest
            end                             % scale. This is used for noise estimation.
 
        end                                 % ... and process the next scale
 
        % Now calculate the phase symmetry measure.
 
        if polarity == 0     % look for 'white' and 'black' spots
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    + abs(real(EO{s,o})) - abs(imag(EO{s,o}));
            end
            
        elseif polarity == 1  % Just look for 'white' spots
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    + real(EO{s,o}) - abs(imag(EO{s,o}));
            end
            
        elseif polarity == -1  % Just look for 'black' spots
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    - real(EO{s,o}) - abs(imag(EO{s,o}));
            end      
             elseif polarity == 2  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    +abs(imag(EO{s,o})) + abs(real(EO{s,o}));
            end     
      
                   elseif polarity == -2  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   +abs(imag(EO{s,o})) - (real(EO{s,o}));
            end  
                               elseif polarity == -3  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   +(imag(EO{s,o})) - (real(EO{s,o}));
            end 
            
                               elseif polarity == 3  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   -(imag(EO{s,o})) - (real(EO{s,o}));
            end 
                               elseif polarity == -4  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   + (real(EO{s,o}))-(imag(EO{s,o})) ;
            end 
                                           elseif polarity == 4  % PhaseAsymmetry
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   + (real(EO{s,o}))+(imag(EO{s,o})) ;
            end 
                                           elseif polarity == -5  % PhaseAsymmetry
           for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                   + abs(imag(EO{s,o}))-abs(real(EO{s,o})) ;
                                                      
 
            end 
        end
        
        % Compensate for noise
        % We estimate the noise power from the energy squared response at the
        % smallest scale.  If the noise is Gaussian the energy squared will
        % have a Chi-squared 2DOF pdf.  We calculate the median energy squared
        % response as this is a robust statistic.  From this we estimate the
        % mean.  The estimate of noise power is obtained by dividing the mean
        % squared energy value by the mean squared filter value
        
        medianE2n = median(reshape(abs(EO{1,o}).^2,1,rows*cols));
        meanE2n = -medianE2n/log(0.5);
        estMeanE2n = [estMeanE2n meanE2n];
        
        noisePower = meanE2n/EM_n;                       % Estimate of noise power.
        
        % Now estimate the total energy^2 due to noise
        % Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))
        
        EstSumAn2 = zero;
        for s = 1:nscale
            EstSumAn2 = EstSumAn2+ifftFilterArray{s}.^2;
        end
        
        EstSumAiAj = zero;
        for si = 1:(nscale-1)
            for sj = (si+1):nscale
                EstSumAiAj = EstSumAiAj + ifftFilterArray{si}.*ifftFilterArray{sj};
            end
        end
        
        EstNoiseEnergy2 = 2*noisePower*sum(sum(EstSumAn2)) + 4*noisePower*sum(sum(EstSumAiAj));
        
        tau = sqrt(EstNoiseEnergy2/2);                % Rayleigh parameter
        EstNoiseEnergy = tau*sqrt(pi/2);              % Expected value of noise energy
        EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );
        
        T =  EstNoiseEnergy + k*EstNoiseEnergySigma;  % Noise threshold
        
        % The estimated noise effect calculated above is only valid for the PC_1
        % measure.  The PC_2 measure does not lend itself readily to the same
        % analysis.  However empirically it seems that the noise effect is
        % overestimated roughly by a factor of 1.7 for the filter parameters
        % used here.
        T = T/1.3;    
 
        % Apply noise threshold 
        Energy_ThisOrient = max(Energy_ThisOrient - T, zero); 
                          
        % Update accumulator matrix for sumAn and totalEnergy
        totalSumAn  = totalSumAn + sumAn_ThisOrient;
        totalEnergy = totalEnergy + Energy_ThisOrient;
        
        % Update orientation matrix by finding image points where the energy in
        % this orientation is greater than in any previous orientation (the
        % change matrix) and then replacing these elements in the orientation
        % matrix with the current orientation number.
        
        if(o == 1),
            maxEnergy = Energy_ThisOrient;
        else
            change = Energy_ThisOrient > maxEnergy;
            orientation = (o - 1).*change + orientation.*(~change);
            maxEnergy = max(maxEnergy, Energy_ThisOrient);
        end
        
    end  % For each orientation
    %fprintf('                                   \r');
    
    
%    disp('Mean Energy squared values recorded with smallest scale filter at each orientation');
%    disp(estMeanE2n);
    
    % Normalize totalEnergy by the totalSumAn to obtain phase symmetry
    OPS = totalEnergy ./ (totalSumAn + epsilon);
%     figure; sc(phaseSym); colormap gray
 

[a b] =size (spread{o});
% size (spread{1})

for i = 1:3
    spread{i}=circshift(spread{i},[a/2 b/2]);
end
mm=spread{1};
for i = 1:2
    mm=mm+spread{i+1};
end

%mm= normalize(normalize(OPS)+normalize(im)); figure, imagesc(mm), figure
%figure,imagesc(spread{1}), figure,imagesc(spread{2}),figure,imagesc(spread{3}),figure,imagesc(mm),figure
%figure,imagesc(ds), figure,imagesc(dc),figure,imagesc(dtheta),figure