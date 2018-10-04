function mYy = mask_Yy(Yy) %#codegen

[~, ncol] = size(Yy);
mYy = Yy;
for icol = 1:ncol
    if any(mYy(:,icol)==1)
        mYy(1:find(Yy(:,icol)==1,1),icol)=1;
    else
         mYy(:,icol)=1;
    end
end
