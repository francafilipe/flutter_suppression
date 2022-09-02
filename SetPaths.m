folders = {'aerodynamics'
           'lqr-optimization'
           'data-results'
           'analysis-results'
           'control'
           'model'
            };
		
for i=1:length(folders)
	f = folders{i};
		cd (f)
		addpath(genpath(pwd))
		cd ..
end

clear f folders i
