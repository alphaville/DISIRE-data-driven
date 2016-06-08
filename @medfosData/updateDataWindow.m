function newWin = updateDataWindow(obj,W, newY, newU, output)
%updates data window with predicted values of Y and measured values of U
szy = size(newY);

if szy(1) > szy(2)
    newY = newY';
end

%W.Y = circshift( W.Y, [1 0] );
indices  = obj.yCross(output,:) > 0;
W.Y = circshift( W.Y, [1 0] );
W.Y(1,:) = newY(indices);

indices  = obj.uCross(output,:) > 0;
W.U = circshift( W.U, [1 0] );
W.U(1, :) = newU(indices);


newWin = W;
return