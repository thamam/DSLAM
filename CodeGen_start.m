
%% Mex mexcostfunc
n=3;
T = 10;
ARGS{1} = coder.typeof(0,[inf,10],[1 1]); %XTstar
ARGS{2} = coder.typeof(0,[3,inf],[0 1]); %XLstar
ARGS{3} = coder.typeof(0,[inf,10],[1 1]);%tildPveccmat
ARGS{4} = coder.typeof(0,[inf,inf,10],[1 1 1]);%Zout 
ARGS{5} = coder.typeof(0,[inf,inf,10],[1 1 1]);%ULout
ARGS{6} = coder.typeof(0,[inf,inf,10],[1 1 1]);%URout
ARGS{7} = coder.typeof(0,1,0); %T
ARGS{8} = coder.typeof(0,1,0);%n
ARGS{9} = coder.typeof(0,1,0);%buffersize



codegen mexcostfunc.m -args ARGS 
% coder.typeof( XTstar,[10,inf],[1 1])
%%coder.typeof(example_value, size_vector [100, inf], variable_dims 0/1 for each dim [ 0 1])

%% Mex mexslamgrad
tf = 10;
BufferSize = 500;
ARGS2{1} = coder.typeof(0,[inf,10],[1 1]); %XTstar
ARGS2{2} = coder.typeof(0,[3,inf],[0 1]); %XLstar
ARGS2{3} = coder.typeof(0,[inf,10],[1 1]);%tildPveccmat
ARGS2{4} = coder.typeof(0,[inf,inf,10],[1 1 1]);%Zout 
ARGS2{5} = coder.typeof(0,[inf,inf,10],[1 1 1]);%ULout
ARGS2{6} = coder.typeof(0,[inf,inf,10],[1 1 1]);%URout
ARGS2{7} = coder.typeof(0,1,0); %T
ARGS2{8} = coder.typeof(0,1,0);%buffersize
codegen mexslamgrad.m -args ARGS2
%% Mex Proj2SE3
ARGS3{1} = coder.typeof(0,[inf,inf,10],[1 1 1]);%XT(just pose) 
ARGS3{2} = coder.typeof(0,1,0); %n
ARGS3{3} = coder.typeof(0,1,0); %tf
ARGS3{4} = coder.typeof(0,1,0); %Bufsize

codegen mexproj2SE3_strct.m -args ARGS3

