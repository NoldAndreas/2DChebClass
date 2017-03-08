

function w = weight_FMT_Roselfeld_x(ptsPolLoc)
    global R
	w = -1/R*ptsPolLoc.y1_kv.*cos(ptsPolLoc.y2_kv);
end
