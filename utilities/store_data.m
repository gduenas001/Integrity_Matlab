function store_data(x, P, xtrue)
% add current data to offline storage
global DATA
CHUNK= 5000;
len= size(DATA.path,2);
if DATA.i == len % grow array exponentially to amortise reallocation
    if len < CHUNK, len= CHUNK; end
    DATA.path= [DATA.path zeros(3,len)];
    DATA.true= [DATA.true zeros(3,len)];
%     pack
end
i= DATA.i + 1;
DATA.i= i;
DATA.path(:,i)= x(1:3);
DATA.true(:,i)= xtrue;
DATA.state(i).x= x;
%DATA.state(i).P= P;
DATA.state(i).P= diag(P);

%
%
