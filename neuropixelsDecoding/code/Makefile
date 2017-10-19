CODEDIR=C:\Users\Peter\Documents\MATLAB\WidefieldDecoding\code
MATRUN=matlab -wait -nodisplay -nosplash -nodesktop -r "cd '$(CODEDIR)', <script>, exit"


preproc/%.mat: code/getAllSessions.m
	$(subst <script>,getAllSessions,$(MATRUN))
	
	