library(ggplot2)
library(maps)
library(mapdata)

state <- map_data("state")

idaho <- subset(state, region=="idaho")
counties <- map_data("county")
idaho_county <- subset(counties, region=="idaho")

ca_map <- ggplot(data=idaho, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill=NA) + 
  geom_polygon(data=idaho_county, fill=NA, color="gray80") + 
  geom_polygon(color="black", fill=NA) + 
  xlim(-120, -110)+
  ylab("Latitude") +
  xlab("Longitude") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme_void()


sites <- read.csv('../data/legacy_data_nw_man_screened.csv')

#need to remove duplicates for labels
sites <- sites[,-2]

#remove lower dry creek and hayspur
sites <- sites[c(-13, -20),]

#remove duplicates
sites <- sites[!duplicated(sites), ]

#keithley not removed??
sites <- sites[c(-5),]

#for map purposes dropping one Mann Cr
sites <- sites[c(-11),]

map <- ca_map +   geom_point(data = sites, aes(x = Longitude, y = Latitude, group = Population)) +
  geom_text_repel(data = sites, aes(x = Longitude, y = Latitude, group = Population, label = Population), hjust=1, vjust=0.5)

ggsave('../outputs/figures/map.jpeg', plot = map, width = 20, height = 20, units = "cm")
