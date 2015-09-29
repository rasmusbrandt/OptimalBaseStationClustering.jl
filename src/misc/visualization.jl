function display_clustering(network, partition, title)
    I = get_num_BSs(network); K = get_num_MSs(network)

    fig = PyPlot.figure()
    ax = fig[:add_axes]((0.15,0.15,0.95-0.15,0.95-0.15))

    # BSs
    local line_BS
    for i = 1:I
        pos = network.BSs[i].position
        line_BS = ax[:plot](pos.x, pos.y; marker="o", color="k", markersize=4, linewidth=0, label="BS")
    end

    # MSs
    local line_MS
    for k = 1:K
        pos = network.MSs[k].position
        line_MS = ax[:plot](pos.x, pos.y; marker="o", color="r", markersize=2, linewidth=0, label="MS")
    end

    # Annotations
    ax[:set_xlabel]("x coordinate [m]")
    ax[:set_ylabel]("y coordinate [m]")
    legend = ax[:legend](handles=[line_BS[1], line_MS[1]], loc="upper left", numpoints=1)
    legend_frame = legend[:get_frame]()
    PyPlot.setp(legend_frame, linewidth=0.5)

    # Square plot
    ax[:set_aspect]("equal", "datalim")

    # Show cluster lines
    for block in partition.blocks
        for i in block.elements
            for j in block.elements
                if i == j; continue; end
                pos1 = network.BSs[i].position; pos2 = network.BSs[j].position
                ax[:plot]([pos1.x, pos2.x], [pos1.y, pos2.y]; marker="o", color="DarkGray", markersize=2)
            end
        end
    end

    fig[:suptitle](title)
    display(fig)
end
