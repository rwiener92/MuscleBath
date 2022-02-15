%%% Ideas for MB-MTC interval detection %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find where force transitions neg->pos
ji = [];
for j = 1:length(Ringdat)-1
	if Ringdat(j,2)<0 && Ringdat(j+1,2)>0
		ji = [ji, j];
	end
end
clear j



%% Apply neg->pos loop to finding next motor step within tolerance
step_size = 0.1; % [mm]
mot_res = 0.07; % [mm]
%motor_step_detect = step_size+mot_res;

j_st=[];
for i_st = 1:length(Ringdat)-1
	app_eq = @(Ringdat(i_st+1,1), Ringdat(i_st,1)+step_size, mot_res) abs(Ringdat(i_st+1,1) - (Ringdat(i_st,1)+step_size)) <mot_res;
	Out = app_eq(Ringdat(i_st+1,1), Ringdat(i_st,1)+step_size, mot_res)

	if Out == 1
		j_st = [j_st, i_st];
	end
end
%%%%%%
 

x = 1;
y = 1.02;
tol = 0.05;
app_eq = @(x,y,tol) abs(x-y)<tol;
Out = app_eq(x,y,tol)
>>Out = 1

x = Ringdat(i_st+1,1);
y = Ringdat(i_st,1);


%%%%%
step_size = 0.1; % [mm]
mot_res = 0.07; % [mm]
%motor_step_detect = step_size+mot_res;

app_eq = @(Ringdat(i_st+1,1), Ringdat(i_st,1)+step_size, mot_res) abs(Ringdat(i_st+1,1) - (Ringdat(i_st,1)+step_size)) <mot_res;
Out = app_eq(Ringdat(i_st+1,1), Ringdat(i_st,1)+step_size, mot_res)
Out



