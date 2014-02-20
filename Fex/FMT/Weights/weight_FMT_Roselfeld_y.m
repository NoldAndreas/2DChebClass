function w = weight_FMT_Roselfeld_y(ptsPolLoc)
    global R
	w = -1/R*ptsPolLoc.y1_kv.*sin(ptsPolLoc.y2_kv);
end