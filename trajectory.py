def drawTrajectory(plt, geodetic, plateCarree, trajectory) -> None:
    l = len(trajectory)
    for i in range(l-1):
        s_lon, s_lat = plateCarree.transform_point(trajectory[i][0], trajectory[i][1], geodetic)
        e_lon, e_lat = plateCarree.transform_point(trajectory[i+1][0], trajectory[i+1][1], geodetic)
        plt.plot([s_lon, e_lon], [s_lat, e_lat], color='black', linewidth=1, marker='o', markersize=3, transform=plateCarree)
