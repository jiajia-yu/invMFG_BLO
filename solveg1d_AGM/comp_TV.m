function tvg = comp_TV(g)

tvg = sum(abs(diff(g)));

end