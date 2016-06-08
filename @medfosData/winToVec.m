function T = winToVec( obj, W, indY )
%a matvec operation on window structure
vector = [];

%do a logic opeartion on inputs and sum
%true (1) to get the number of inputs
nReg = obj.yCross(indY, obj.yCross(indY,:) > 0 );
for i = 1:sum( obj.yCross(indY, :) > 0 )
        vector = [vector; W.Y(1:nReg(i),i) ];
end

%outputs
nReg = obj.uCross(indY, obj.uCross(indY,:) > 0 );
for i = 1:sum( obj.uCross(indY, :) > 0 )
        vector = [vector; W.U(1:nReg(i),i) ];
end

T = vector; 
end