function value = checkenv(envVarName,default)
% CHECKENV.m
%
% Takes an environment variable name and sets that, if it's there
%
% USAGE
% seed = checkenv('SEEDVALUE',sum(clock*100));
%  - this checks the environment for SEEDVALUE, and sets seed
%    to this env value if its there. sets to sum(...) if it's not
%
% NOTE
% only useful for float/int and string

% set default
value = default;

% grab env variable
tmp = getenv(envVarName);

% if its not empty, set!
if ~isempty(tmp)
	if ischar(default)
		value = tmp;
	else
		value = str2num(tmp);
	end
end

end %function