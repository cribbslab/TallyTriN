

class style(object):

    def __init__(self, ):
        pass

    def deephelicon(self, boxplot_handles, palette=[]):
        """

        ..  Description:
            ------------
            1. bplot_handle offers six properties to adjust, including:
            1) 'whiskers', 2) 'caps', 3) 'medians', 4) 'fliers', 5) 'means', 6) 'boxes'.

            2. style scheme used in Figure S1 in the DeepHelicon (link 1) paper.

        ..  See:
            ----
            1. DeepHelicon: https://doi.org/10.1016/j.jsb.2020.107574.

        :param boxplot_handles: dict - plt boxplot obs
        :param palette: 1d list
        :return:
        """
        for bplot_handle, color in zip(boxplot_handles, palette):
            for patch in bplot_handle['boxes']:
                patch.set(color='black', linewidth=1.1)
                patch.set(facecolor='white', alpha=0.3)
            for whisker in bplot_handle['whiskers']:
                whisker.set(color='black', linewidth=1)
            for cap in bplot_handle['caps']:
                cap.set(color='black', linewidth=1)
            for median in bplot_handle['medians']:
                median.set(color='black', linewidth=1)
            for flier in bplot_handle['fliers']:
                flier.set(marker='o', color='y', alpha=0.5)
        return

    def plain(self, boxplot_handles, palette, edge_width=1):
        """

        ..  Description:
            ------------
            bplot_handle offers six properties to adjust, including:
            1) 'whiskers', 2) 'caps', 3) 'medians', 4) 'fliers', 5) 'means', 6) 'boxes'.


        :param boxplot_handles: dict - plt boxplot obs
        :param palette: 1d list
        :param edge_width:
        :return:
        """
        for bplot_handle, color in zip(boxplot_handles, palette):
            for patch in bplot_handle['boxes']:
                # b.set_alpha(0.6)
                # patch.set_edgecolor(color, alpha=0.3)
                patch.set(
                    edgecolor=color,
                    alpha=0.7,
                )
                # patch.set_facecolor('white')
                patch.set_linewidth(edge_width)
                patch.set(
                    facecolor='white'
                )
            for whisker in bplot_handle['whiskers']:
                whisker.set(
                    color=color,
                    linewidth=2,
                    linestyle='dashed', # the same as '--'
                    alpha=0.5,
                )
            for cap in bplot_handle['caps']:
                cap.set(color=color, linewidth=1)
            for median in bplot_handle['medians']:
                median.set(color=color, linewidth=1)
            for flier in bplot_handle['fliers']:
                flier.set(marker='o', color='y', alpha=0.5)

    def plaincolor(self, plt, boxplot_handles, palette=[]):
        """
        ..  Description:
            ------------
            only adjust uniform colors each boxplot.

        :param plt:
        :param boxplot_handles: dict - plt boxplot obs
        :param palette: 1d list
        :return:
        """
        def set(bp_handle, color):
            plt.setp(bp_handle['boxes'], color=color)
            plt.setp(bp_handle['whiskers'], color=color)
            plt.setp(bp_handle['caps'], color=color)
            plt.setp(bp_handle['medians'], color=color)
        for i_bp, c in zip(boxplot_handles, palette):
            set(i_bp, c)
        return