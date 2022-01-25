function [ caCode ] = gen_ca_code( filePath, prnList )
%GEN_CA_CODE Returns C/A code for specified PRN list.
%   Returns a MxN matrix where each row corresponds to
% a PRN where M is length of input argument list and N
% is equal to the number of chips (1023).
fileName = [filePath, 'caCode.mat'];
genCaCode = true;

if exist(fileName, 'file') == 2
    caCode=load(fileName);
    if (length(prnList) == 1 && any(ismember(prnList, caCode.prnList))) || ...
            all(caCode.prnList == prnList)
        caCode=caCode.caCode(prnList, :);
        genCaCode = false;
    else
        genCaCode = true; 
    end
end

if genCaCode
        
    %--- Make the code shift array. The shift depends on the PRN number -------
    % The g2s vector holds the appropriate shift of the g2 code to generate
    % the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
    g2s = [  5,   6,   7,   8,  17,  18, 139, 140, 141, 251, ...
        252, 254, 255, 256, 257, 258, 469, 470, 471, 472, ...
        473, 474, 509, 512, 513, 514, 515, 516, 859, 860, ...
        861, 862 ... end of shifts for GPS satellites
        ... Shifts for the ground GPS transmitter are not included
        ... Shifts for EGNOS and WAAS satellites (true_PRN = PRN + 87)
        145, 175,  52,  21, 237, 235, 886, 657, ...
        634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, ...
        386];
    
    
    % Parameters
    caCodeLength = 1023;
    numPrns = length(prnList);
    
    % C/A code buffer
    caCode = zeros(numPrns, caCodeLength);
    
    % Iterate over PRN list
    for prnIdx = 1:numPrns
        prn = prnList(prnIdx);
        % Generate all G1 & G2 signal chips
        g1 = zeros(1, caCodeLength);
        g2 = zeros(1, caCodeLength);
        reg1 = -1*ones(1, 10);
        reg2 = -1*ones(1, 10);
        for i=1:caCodeLength
            g1(i)     = reg1(10);
            g2(i)     = reg2(10);
            saveBitG1 = reg1(3)*reg1(10);
            saveBitG2 = reg2(2)*reg2(3)*reg2(6)*reg2(8)*reg2(9)*reg2(10);
            reg1      = circshift(reg1, 1);
            reg2      = circshift(reg2, 1);
            reg1(1)    = saveBitG1;
            reg2(1)    = saveBitG2;
        end
        
        % Form single sample C/A code by multiplying G1 and G2
        g2shift = g2s(prn);
        caCode(prnIdx, :) = -(g1 .* circshift(g2, g2shift));
    end
    save(fileName, 'caCode', 'prnList');
    
end
end




