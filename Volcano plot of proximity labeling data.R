rm(list=ls())
library(ggplot2)

######### load data ############################################################################################
counter.receptor=read.csv('Counter receptors', row.names = 1)

######### Volcano plot GOCC #####################################################################################

# all black
volcano.plot=ggplot(counter.receptor,aes(LFQ.diff,-log10(LFQ.p)))+
  geom_point(size=0.5) + 
  xlab('log2 difference') + 
  ylab('-log10(p value)') +
  scale_x_continuous(breaks = seq(-10,10,by=1),limits = c(-4.5,4.5)) +
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  labs(title = 'SnAct All Proteins Black') + 
  theme(text = element_text(size=10), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5, size=10)) +
  geom_vline(xintercept = c(-1,1),linetype='dotted', color='blue') +  geom_hline(yintercept = -log10(0.05),linetype='dotted', color='blue')

ggsave('Volcano all black', device = 'tiff', dpi='retina',width = 5, height = 3, units = 'in')


# nucleus red
volcano.plot+
geom_point(data=subset.data.frame(counter.receptor, GOCC.annotation=='nucleus'), 
           aes(LFQ.diff,-log10(LFQ.p)),color='red', size=0.5) + 
           labs(title = 'SnAct Nucleus Red')

ggsave('Volcano nucleus red', device = 'tiff', dpi='retina',width = 5, height = 3, units = 'in')

# cytoplasm red
volcano.plot+
  geom_point(data=subset.data.frame(counter.receptor, GOCC.annotation=='cytoplasm'), 
             aes(LFQ.diff,-log10(LFQ.p)),color='red', size=0.5) + 
  labs(title = 'SnAct Cytoplasm Red')

ggsave('Volcano cytoplasm red', device = 'tiff', dpi='retina',width = 5, height = 3, units = 'in')

# other red
volcano.plot+
  geom_point(data=subset.data.frame(counter.receptor, GOCC.annotation=='other'), 
             aes(LFQ.diff,-log10(LFQ.p)),color='red', size=0.5) + 
  labs(title = 'SnAct Other Red')

ggsave('Volcano other red', device = 'tiff', dpi='retina',width = 5, height = 3, units = 'in')

# plasma membrane red
volcano.plot+
  geom_point(data=subset.data.frame(counter.receptor, GOCC.annotation=='plasma membrane'), 
             aes(LFQ.diff,-log10(LFQ.p)),color='red', size=0.5) + 
  labs(title = 'SnAct Plasma Membrane Red')

ggsave('Volcano plasma membrane red', device = 'tiff', dpi='retina',width = 5, height = 3, units = 'in')