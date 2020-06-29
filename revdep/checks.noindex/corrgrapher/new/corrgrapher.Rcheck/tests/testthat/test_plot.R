context('plot working properly')

df <- as.data.frame(datasets::Seatbelts)
df <- df[,1:(ncol(df)-1)]

expect_is(plot(corrgrapher(df)), 'visNetwork')

