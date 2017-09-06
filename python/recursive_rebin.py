# def recursive_rebin(x, y, y_err):
#     # def grouper(iterable, n, fillvalue=None):
#     #     args = [iter(iterable)] * n
#     #     return izip_longest(*args, fillvalue=fillvalue)

#     # Discard points with zero error
#     print x, len(x)
#     _x, _y, _y_err = zip*([_x, _y, _y_err for _x, _y, _yerr in zip(x, y, y_err) if y_err != 0])
#     print "->", x, len(x)


#     # # if len(x) % rebin != 0:
#     # #     print "Rebinning {} points by factor {} will leave ungrouped points!".format(len(x), rebin)

#     # # while True:
#     # fractional_errors = [ _y_err / _y for _y, _y_err in zip(x, y, y_err) if y_err != 0]
#     # idx_highest_error = datapoints.index(max(datapoints))


#     # for group in grouper(zip(x, y, y_err), rebin):
#     #     try:
#     #         group = [x for x in group if x is not None]
#     #         _reciprocal_sqrd_err = sum([1.0/ (elem[2]**2) for elem in group])
#     #         _x.append(sum([elem[0] for elem in group]) / len(group))
#     #         _y.append(sum([elem[1] / (elem[2]**2) for elem in group]) / _reciprocal_sqrd_err)
#     #         _y_err.append(1.0 / np.sqrt(_reciprocal_sqrd_err))
#     #     except:
#     #         print group
#     #         raise
#     # return _x, _y, _y_err







# # def combine_bins(bin1, bin2):
# #     x = (bin1[0] + bin2[0]) / 2
# #     # y = (bin1.y/bin1.y_err + bin2.y/bin2.y_err) / (1.0 / bin1.y_err + 1.0 / bin2.y_err)
# #     y = (bin1.y/bin1.y_err + bin2.y/bin2.y_err) / (1.0 / bin1.y_err + 1.0 / bin2.y_err)

# # # if self.value()[1] == 0 : self.y_err = 1e-99
# # #     if other.value()[1] == 0 : other.y_err = 1e-99
# # #     x_low, x_high = min( self.range()[0], other.range()[0] ), max( self.range()[1], other.range()[1] )
# # #     y = ( self.value()[0]/self.value()[1] + other.value()[0]/other.value()[1] ) / ( 1.0/self.value()[1] + 1.0/other.value()[1] )
# # #     y_err = math.sqrt( 2.0 / ( 1.0/self.value()[1] + 1.0/other.value()[1] )**2 )

# # def
# # def recursive_rebin(dataset):
# #     bins = [(_x, _y, _y_err, abs(_y / _y_err)) for _x, _y, _y_err in dataset]

# #     # while True:
# #     #     # Get lowest significance bin and exit if this is above the cut
# #     #     significance_ordered_bins = sorted(bins, key=lambda x: x[3])
# #     #     lowest_significance, idx_lowest_significance = significance_ordered_bins[0][3], significance_ordered_bins[0][1]
# #     #     if lowest_significance > significance_cut : break

# #         # # Get bins to merge
# #         # if len( decorated_bins ) == 1 : break
# #         # adjacent_bins = []
# #         # if idx_lowest_significance > 0 : adjacent_bins.append( decorated_bins[idx_lowest_significance-1] )
# #         # if idx_lowest_significance < len(decorated_bins)-1 : adjacent_bins.append( decorated_bins[idx_lowest_significance+1] )
# #         # idx_merge_bin = sorted( adjacent_bins, key=operator.itemgetter(2) )[0][1]

# # # # def significance(self) :
# # # #     difference, error = abs( 1.0 - self.y ), abs( self.y_err )
# # # #     if difference < 1e-10 or error < 1e-10 : return 1e99
# # # #     return difference / error

# # # def recursive_1D_rebin(x, y, y_err) :
# # # #   input_sets = [ DataPoint( input_hist.GetBinLowEdge(idx_bin), input_hist.GetBinLowEdge(idx_bin+1), input_hist.GetBinContent(idx_bin), input_hist.GetBinError(idx_bin) ) for idx_bin in range(1,input_hist.GetNbinsX()+1) ]

# # # #   # Set up initial bin list
# # # #   decorated_bins = [ [ input_set, idx, input_set.significance() ] for idx, input_set in enumerate(input_sets) ]

# # #     bins = []

# # #     while True :
# # #         # print [ (point.range(),point.value()) for point in [ bin[0] for bin in decorated_bins ] ]
# # #         # Get lowest significance bin and exit if this is above the cut
# # #         significance_ordered_bins = sorted( decorated_bins, key=operator.itemgetter(2) )
# # #         lowest_significance, idx_lowest_significance = significance_ordered_bins[0][2], significance_ordered_bins[0][1]
# # #         if lowest_significance > significance_cut : break

# # # #     # Get bins to merge
# # # #     if len( decorated_bins ) == 1 : break
# # # #     adjacent_bins = []
# # # #     if idx_lowest_significance > 0 : adjacent_bins.append( decorated_bins[idx_lowest_significance-1] )
# # # #     if idx_lowest_significance < len(decorated_bins)-1 : adjacent_bins.append( decorated_bins[idx_lowest_significance+1] )
# # # #     idx_merge_bin = sorted( adjacent_bins, key=operator.itemgetter(2) )[0][1]

# # # #     # Create merged data point and replace the existing bins
# # # #     new_data_point = decorated_bins[idx_lowest_significance][0] + decorated_bins[idx_merge_bin][0]
# # # #     decorated_bins[ idx_lowest_significance ] = [ new_data_point, idx_lowest_significance, new_data_point.significance() ]
# # # #     del decorated_bins[ idx_merge_bin ]
# # # #     for idx, bin in enumerate(decorated_bins) : bin[1] = idx

# # # #   # Get bin edges for output histogram
# # # #   bin_edges, seen = sum( [ [decorated_bin[0].range()[0], decorated_bin[0].range()[1]] for decorated_bin in decorated_bins ], [] ), set()
# # # #   bin_edges = array.array( 'd', [ x for x in bin_edges if x not in seen and not seen.add(x) ] )

# # # #   output_hist = ROOT.TH1D( str(uuid.uuid4()), input_hist.GetTitle(), len(bin_edges)-1, bin_edges )
# # # #   if len(decorated_bins) == 1 and abs(1.0-decorated_bins[idx][0].value()[0]) < decorated_bins[idx][0].value()[1] :
# # # #     decorated_bins[idx][0].y = 1.0 # set content to 1.0 if bin is compatible with 1.0
# # # #   [ output_hist.SetBinContent( idx+1, decorated_bins[idx][0].value()[0] ) for idx in range(len(decorated_bins)) ]
# # # #   [ output_hist.SetBinError( idx+1, decorated_bins[idx][0].value()[1] ) for idx in range(len(decorated_bins)) ]
# # # #   # Set end points
# # # #   output_hist.SetBinContent( 0, input_hist.GetBinContent(0) ); output_hist.SetBinContent( output_hist.GetNbinsX()+1, input_hist.GetBinContent(input_hist.GetNbinsX()+1) );
# # # #   output_hist.SetBinError( 0, input_hist.GetBinError(0) ); output_hist.SetBinError( output_hist.GetNbinsX()+1, input_hist.GetBinError(input_hist.GetNbinsX()+1) )
# # # #   return output_hist