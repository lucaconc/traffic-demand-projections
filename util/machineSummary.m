function s = machineSummary()
%MACHINESUMMARY One-line description of the computer used for timing.
%   s = machineSummary() returns e.g. 'MATLAB 24.2 (R2024b), maca64, 10 cores,
%   Apple M1 Pro'. It is stored with the results, so that computation times
%   are reported for the machine that produced them.
parts = {};
try
    v = ver('MATLAB'); parts{end+1} = sprintf('MATLAB %s %s',v.Version,v.Release);
catch
end
parts{end+1} = computer('arch');
try
    parts{end+1} = sprintf('%d cores',feature('numcores'));
catch
end
cpu = '';
try
    if ispc
        cpu = strtrim(getenv('PROCESSOR_IDENTIFIER'));
    elseif ismac
        [st,out] = system('sysctl -n machdep.cpu.brand_string'); if st==0, cpu = strtrim(out); end
    else
        [st,out] = system('grep -m1 "model name" /proc/cpuinfo | cut -d: -f2'); if st==0, cpu = strtrim(out); end
    end
catch
end
if ~isempty(cpu), parts{end+1} = cpu; end
s = strjoin(parts,', ');
end
