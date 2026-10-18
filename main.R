# Kevin Hon Yin Hau
# k-Bit Recurrence Model (kBRM) -- R Implementation
# Includes: absorbing barrier (maxPay), cyclostationary extension (mth_num),
#           full joint distribution output, and decline curve calculation.

library(dplyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)


## ---- paidOrNotNext ---------------------------------------------------
## Single-step propagation. Takes the current-step distribution table and
## the parameter table for this step, returns the next-step distribution.
## By operating on the current step only (not full history), memory stays
## O(2^k * min(t, T)) rather than growing quadratically with t.
## payCodeTbl is passed separately because it can change at each step
## (cyclostationary case) or remain fixed (static case).

paidOrNotNext <- function(tbl, payCodeTbl, maxPay = 999) {
  
  ## Extract the most recent step and join payment probabilities.
  ## Shift the state code left by one position (drop oldest bit).
  tbl_last <- tbl %>%
    filter(numStep == max(tbl$numStep)) %>%
    inner_join(payCodeTbl, c("code" = "code")) %>%
    mutate(code = substr(code, 2, nchar(tbl$code[1])))
  
  out <- tbl %>%
    union_all(
      ## Paid branch: append '1', increment payNum (capped at maxPay),
      ##              multiply probability by p(s).
      tbl_last %>%
        mutate(code   = paste(code, "1", sep = ""),
               payNum = if_else(payNum >= maxPay, maxPay, payNum + 1L),
               prob   = prob * paidPerc) %>%
        union_all(
          ## Not-paid branch: append '0', payNum unchanged,
          ##                  multiply probability by q(s) = 1 - p(s).
          tbl_last %>%
            mutate(code = paste(code, "0", sep = ""),
                   prob = prob * notPaidPerc)
        ) %>%
        mutate(numStep = numStep + 1L) %>%
        select(-paidPerc, -notPaidPerc)
    )
  
  ## Collapse: sum probabilities over identical (numStep, code, payNum) tuples.
  ## This is the key operation that keeps the state space O(2^k) rather
  ## than O(2^t) -- many paths share the same state and payment count.
  out %>%
    group_by(numStep, code, payNum) %>%
    summarise(prob = sum(prob)) %>%
    ungroup()
}


## ---- paidOrNotPaid ---------------------------------------------------
## Full propagation for one starting state over numMth steps.
## Handles both static (no mth_num column) and cyclostationary
## (mth_num column, cycle length = max(mth_num)) parameter tables.
## Returns the joint distribution (startCode, numStep, payNum, prob)
## for all steps 1 to numMth.

paidOrNotPaid <- function(startCode, payCodeTbl, numMth,
                          start_mth_num = 1, maxPay = 999) {
  
  chk_cycle <- any(names(payCodeTbl) == "mth_num")
  
  ## Pre-select columns for static case (avoids repeated select inside loop).
  if (!chk_cycle)
    pay_code <- payCodeTbl %>% select(code, paidPerc, notPaidPerc)
  
  ## Initialise: step 0, starting state, zero payments, probability = 1.
  tbl <- data.frame(numStep = 0L, code = startCode, payNum = 0L, prob = 1)
  
  for (i in seq_len(numMth)) {
    
    ## Cyclostationary: look up the parameter set for the current calendar
    ## month. The modular index maps any forecast horizon back to the
    ## correct month in the cycle.
    if (chk_cycle) {
      mth_idx  <- ((start_mth_num + i) %% max(payCodeTbl$mth_num)) + 1L
      pay_code <- payCodeTbl %>%
        filter(mth_num == mth_idx) %>%
        select(code, paidPerc, notPaidPerc)
    }
    
    tbl <- paidOrNotNext(tbl, pay_code, maxPay)
  }
  
  ## Remove the step-0 initialisation row and add the startCode label.
  tbl %>%
    filter(numStep > 0) %>%
    mutate(startCode = startCode) %>%
    group_by(startCode, numStep, payNum) %>%
    summarise(prob = sum(prob)) %>%
    ungroup()
}


## ---- getPymtPercTbl --------------------------------------------------
## Runs paidOrNotPaid for every unique starting state in payCodeTbl and
## combines the results. Returns the full joint distribution
## P(N_t = m, S_t = s) for all starting states, steps, and payment counts.
## This table is the foundation for all downstream analysis: expected
## value, variance, confidence intervals, and quantiles are all derived
## from it with no additional computation.

getPymtPercTbl <- function(payCodeTbl, numMth,
                           start_mth_num = 1, maxPay = 999) {
  
  unique_codes <- unique(payCodeTbl$code)
  
  result <- paidOrNotPaid(unique_codes[1], payCodeTbl,
                          numMth, start_mth_num, maxPay)
  
  for (i in seq_along(unique_codes)[-1]) {
    result <- union_all(result,
                        paidOrNotPaid(unique_codes[i], payCodeTbl,
                                      numMth, start_mth_num, maxPay))
  }
  
  result
}




notPaidWithinMth <- function(startCode, payCodeTbl, numMth,
                             start_mth_num = 1, maxPay = 999) {
  
  chk_cycle <- any(names(payCodeTbl) == "mth_num")
  
  ## Pre-select columns for static case (avoids repeated select inside loop).
  if (!chk_cycle)
    pay_code <- payCodeTbl %>% select(code, paidPerc, notPaidPerc)
  
  ## Initialise: step 0, starting state, zero payments, probability = 1.
  tbl <- data.frame(numStep = 0L, code = startCode, payNum = 0L, prob = 1)
  
  for (i in seq_len(numMth)) {
    
    ## Cyclostationary: look up the parameter set for the current calendar
    ## month. The modular index maps any forecast horizon back to the
    ## correct month in the cycle.
    if (chk_cycle) {
      mth_idx  <- ((start_mth_num + i) %% max(payCodeTbl$mth_num)) + 1L
      pay_code <- payCodeTbl %>%
        filter(mth_num == mth_idx) %>%
        select(code, paidPerc, notPaidPerc)
    }
    
    tbl <- paidOrNotNext(tbl, pay_code, maxPay)
  }
  
  ## Remove the step-0 initialisation row and add the startCode label.
  tbl %>%
    filter(numStep > 0) %>%
    mutate(startCode = startCode) %>%
    group_by(startCode, numStep, payNum) %>%
    summarise(prob = sum(prob)) %>%
    ungroup()
}










## ---- getModelDeclineCurve --------------------------------------------
## Post-processing only: derives the expected cumulative payments E[N_t]
## and the marginal expected payment (decline curve) D(t) = E[N_t] - E[N_{t-1}]
## from the getPymtPercTbl() output. No further computation is needed.
## D(t) represents the probability-weighted expected payment in month t.
## Under the absorbing barrier model, D(t) -> 0 as t -> infinity.

getModelDeclineCurve <- function(pymtPercTbl) {
  
  result <- pymtPercTbl %>%
    group_by(startCode, numStep) %>%
    summarise(expectedPayt = sum(payNum * prob)) %>%
    ungroup()
  
  result %>%
    left_join(
      result %>%
        mutate(numStep = numStep + 1L, expectedPaytLast = expectedPayt) %>%
        select(startCode, numStep, expectedPaytLast),
      by = c("startCode", "numStep")
    ) %>%
    mutate(declineCurve = expectedPayt - coalesce(expectedPaytLast, 0)) %>%
    select(-expectedPaytLast)
}


## ---- getCodeTbl ------------------------------------------------------
## Helper: converts a raw count table (code, cnt, output) into the
## parameter table (code, total, paidPerc, notPaidPerc) required by
## getPymtPercTbl(). paidPerc = output / cnt is the empirical MLE of
## P(payment | state).

getCodeTbl <- function(dt) {
  data.frame(
    code        = dt$code,
    total       = dt$cnt,
    paidPerc    = dt$output / dt$cnt,
    notPaidPerc = 1 - dt$output / dt$cnt,
    stringsAsFactors = FALSE
  ) %>% arrange(code)
}


## ---- Example usage ---------------------------------------------------
## k = 2 parameter table. Each row is one history state:
##   code '11' = paid, paid   (p = 0.90)
##   code '10' = paid, missed (p = 0.50)
##   code '01' = missed, paid (p = 0.30)
##   code '00' = missed, missed (p = 0.05)

x <- data.frame(
  code   = c("00", "01", "10", "11"),
  cnt    = c(100,  500,  500,  1000),
  output = c( 30,  300,  200,  900)
)
y <- getCodeTbl(x)

example <- getPymtPercTbl(y, 24)

example_result <- getModelDeclineCurve(example)







library(ggplot2)

# =====================================================================
# PLOTS USED IN THE PAPER
# Fig. 1 (shift-register example) and Fig. 2 (de Bruijn diagram) are not
# drawn in R (see the PNG -> EPS conversion at the end of this file).
# Figs. 3-6 are produced below, in the order Fig. 5, 4, 3, 6.
# =====================================================================

# ---- Paper Fig. 5 (fig5.eps), Section 4.3 (label fig:curves) --------
# Payment curves D(t) for each starting state {00, 01, 10, 11} (k = 2,
# 24 months): a decline curve from '11' and a growth curve from '00'.
ggplot(example_result, 
       aes(x = numStep, y = declineCurve),
       breaks = c(5, 10, 15, 20)) + geom_line() + 
  scale_y_continuous(breaks = (1:5)/5, labels = scales::percent_format(), limits = c(0.3, 0.95)) + 
  labs(y = "Expected Payment", x = "Num of Months", 
       title = "") +
  theme_bw() +
  facet_wrap(. ~ startCode, nrow=2) + 
  theme(legend.position="none", 
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10)) 

#------------------------------------------------------------------------------------------

x2 <- x %>% mutate(proportion = cnt / sum(x$cnt))

example_all <-   example %>% 
  inner_join(x2 , by = c("startCode" = "code")) %>% 
  group_by(numStep, payNum) %>% summarise(prob = sum(prob * proportion))

# ---- NOT USED: draft PMF plots ---------------------------------------
# The plots below (PMF by starting state at t = 6, 12, 18, 24, and a
# single marginal PMF at t = 24) were exploratory drafts. They are
# superseded by the 4-panel Fig. 4 further down, so they are commented out.
# Draft: PMF by starting state, t = 6
# ggplot(example %>% filter(numStep == 6), 
#        aes(x = payNum, y = prob, group = startCode)) + 
#   geom_bar(stat="identity") + facet_wrap(. ~ startCode) +
#   labs(y = "p.m.f.", x = "Num of Payments at t = 6", 
#        title = "") +
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.6)) + 
#   theme_bw()+ 
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10)) 



# Draft: PMF by starting state, t = 12
# ggplot(example %>% filter(numStep == 12), 
#        aes(x = payNum, y = prob, group = startCode)) + 
#   geom_bar(stat="identity") + facet_wrap(. ~ startCode) +
#   labs(y = "p.m.f.", x = "Num of Payments at t = 12", 
#        title = "") +
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.3)) + 
#   theme_bw()+ 
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10)) 


# Draft: PMF by starting state, t = 18
# ggplot(example %>% filter(numStep == 18), 
#        aes(x = payNum, y = prob, group = startCode)) + 
#   geom_bar(stat="identity") + facet_wrap(. ~ startCode) +
#   labs(y = "p.m.f.", x = "Num of Payments at t = 18", 
#        title = "") +
#   scale_y_continuous(labels = scales::percent_format()) + 
#   theme_bw()+ 
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10)) 

# Draft: PMF by starting state, t = 24
# ggplot(example %>% filter(numStep == 24), 
#        aes(x = payNum, y = prob, group = startCode)) + 
#   geom_bar(stat="identity") + facet_wrap(. ~ startCode) +
#   labs(y = "p.m.f.", x = "Num of Payments at t = 24", 
#        title = "") +
#   scale_y_continuous(labels = scales::percent_format()) + 
#   theme_bw()+ 
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10)) 
# Draft: marginal PMF at t = 24 only
# ggplot(example_all %>% filter(numStep == 24),
#        aes(x = payNum, y = prob)) +
#   geom_bar(stat="identity", fill="#aaaaaa") +
#   labs(y = "p.m.f.", x = "Num of Payments at t = 24",
#        title = "") +
#   scale_y_continuous(labels = scales::percent_format()) +
#   theme_bw()+
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10))

# ---- Paper Fig. 4 (fig4.eps), Section 3.3 (label fig:pmf) -----------
# PMF P(N_t = m) of the cumulative payment count at t = 6, 12, 18, 24,
# marginalised over starting states by portfolio proportion (example_all).
ggplot(example_all %>% filter(numStep %in% c(6,12,18,24)) %>% 
         mutate(numStep = factor(paste0("t = ", numStep), 
                                 levels = paste0("t = ", (1:4)*6))), 
       aes(x = payNum, y = prob, group = numStep)) + 
  geom_bar(stat="identity", fill="#aaaaaa") + 
  facet_wrap(. ~ numStep, scales = "free") +
  labs(y = "p.m.f.", x = "Num of Payments at t", 
       title = "") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_bw()+ 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10)) 


#--------------------------------------------------------









y2 <- y %>% mutate(mth_num = 1) %>% union(
  y %>% mutate(paidPerc = paidPerc *0.9, notPaidPerc = 1-paidPerc, mth_num = 2)
)

## Step 1: compute full joint distribution over 5 months (no maturity limit).
## Step 2: derive decline curve from the distribution table.
print(data.frame(getModelDeclineCurve(getPymtPercTbl(y, 5))))

print(data.frame(getModelDeclineCurve(getPymtPercTbl(y2, 5))))



## Confidence intervals -- no additional modelling needed.
## Marginalise over states, compute cumulative probabilities, read quantiles.
pymtDist <- getPymtPercTbl(y, numMth = 10, maxPay = 4)

marginal <- pymtDist %>%
  group_by(startCode, numStep, payNum) %>%
  summarise(prob = sum(prob)) %>%
  group_by(startCode, numStep) %>%
  arrange(payNum) %>%
  mutate(cum_prob = cumsum(prob))

## 95% confidence interval: Q_0.025 and Q_0.975
ci <- marginal %>%
  summarise(
    lower_95 = min(payNum[cum_prob >= 0.025]),
    upper_95 = min(payNum[cum_prob >= 0.975])
  )
print(ci)











View(example %>% arrange(numStep))

#----------------------------------------------------------

# ---- Paper Fig. 3 (fig3.eps), Section 3.2 (label fig:complexity) ----
# Terms propagated at each step t: kBRM, at most 2^k (t+1) (flat),
# versus explicit path enumeration, 2^k 2^t (exponential). Log scale.
ggplot(example %>% group_by(numStep) %>% 
         summarise(total = n()) %>% ungroup() %>% 
         mutate(type = "kBRM (flat)") %>% 
         union(data.frame(numStep=1:24, total = (2^2)*2^(1:24), 
                          type="naive binary tree (exponential)")),
       aes(numStep, total, colors = "black", linetype = type)) + geom_line(linewidth = 0.8) +
  geom_text(data = NULL, aes(x = 15, y = 200, label = "italic(2^k*(t+1))"), 
            size = 4, family = "serif", parse=TRUE, colour = "gray30") + 
  geom_text(data = NULL, aes(x = 15, y = 40000, label = "italic(2^k*2^t)"), 
            size = 4, family = "serif", parse=TRUE, colour = "gray30") + 
  scale_y_log10() +
  labs(y="Total outputs (k = 2)", x = "Num of steps (t)") +
  theme_bw() + 
  theme(legend.position="bottom", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank()) 


#-----------------------------------------------------------------------
# ---- Paper Fig. 6 (fig6.eps), Section 5.3 (label fig:cyclo) ---------
# Payment curve D(t): static kBRM vs cyclostationary kBRM with a 12-month
# seasonal cycle, same k = 2 base parameters.


y12 <- y %>% 
  cross_join(data.frame(mth_num=1:12)) %>% 
  inner_join(
    data.frame(mth_num = 1:12, 
               cycle_adj = 0.95 + 0.05 *sin(2*(0:11)*pi/12))
  ) %>% mutate(paidPerc = paidPerc * cycle_adj, 
               notPaidPerc = 1 - paidPerc) %>% 
  arrange(mth_num, code) 


y12 <- y %>% mutate(mth_num = 1) %>% 
  union(y %>% mutate(mth_num = 2)) %>% 
  union(y %>% mutate(mth_num = 3)) %>% 
  union(y %>% mutate(mth_num = 4)) %>% 
  union(y %>% mutate(mth_num = 5)) %>% 
  union(y %>% mutate(mth_num = 6)) %>% 
  union(y %>% mutate(mth_num = 7)) %>% 
  union(y %>% mutate(mth_num = 8)) %>% 
  union(y %>% mutate(mth_num = 9)) %>% 
  union(y %>% mutate(mth_num = 10)) %>% 
  union(y %>% mutate(mth_num = 11)) %>% 
  union(y %>% mutate(mth_num = 12)) %>% 
  inner_join(
    data.frame(mth_num = 1:12, 
               cycle_adj = 0.95 + 0.05 *sin(2*(0:11)*pi/12))
  ) %>% mutate(paidPerc = paidPerc * cycle_adj, 
               notPaidPerc = 1 - paidPerc)

example_cycle12 <- getPymtPercTbl(y12, 24)

example_result12 <- getModelDeclineCurve(example_cycle12)

compareResult <- example_result %>% 
  inner_join(x2, by=c("startCode"="code")) %>% 
  group_by(numStep) %>% 
  summarise(declineCurve = sum(declineCurve * proportion)) %>% 
  mutate(type = "Original kBRM (k=2)") %>% 
  union(
    example_result12 %>% 
      inner_join(x2, by=c("startCode"="code")) %>% 
      group_by(numStep) %>% 
      summarise(declineCurve = sum(declineCurve * proportion)) %>% 
      mutate(type = "Cyclostationary (k=2 with 12-months)")
  ) %>% mutate(
    type = factor(type, 
                  levels = c("Original kBRM (k=2)", "Cyclostationary (k=2 with 12-months)"))
  )

ggplot(compareResult, aes(x=numStep, y=declineCurve, colors = "black", linetype = type)) + 
  geom_line(, linewidth = 0.8) + labs(y="p.m.f", x = "Num of steps (t)") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.5,0.8)) + 
  theme_bw() + 
  theme(legend.position="bottom", 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank()) 

# -----------------------------------------------------------------------------------



# ---- PNG -> EPS conversion (hand-made diagrams, not R plots) ---------
# These appear to be the sources of paper Fig. 1 (shift-register example)
# and Fig. 2 (de Bruijn diagram). Check which of 'fig 0' / 'fig 1' maps
# to fig1.eps / fig2.eps.
library(magick)
fig0 <- image_read("C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 0.png")
image_write(fig0, path = "C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 0.eps", format = "eps")

fig1 <- image_read("C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 1.png")
image_write(fig1, path = "C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 1.eps", format = "eps")

# NOT USED: old test conversion of fig8.png. The final Fig. 3 is the
# ggplot above (fig3.eps).
# fig3 <- image_read("C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig8.png")
# image_write(fig3, path = "C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 3_test.eps", format = "eps")


img <- image_read("C:\\Users\\HYH\\Documents\\HON\\work\\R\\Dynamic Stochastic\\fig\\fig 0.eps") 
plot(img)
