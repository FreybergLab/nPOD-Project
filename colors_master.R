# =============================================================================
# colors_master.R
# Master color palette for all figures — source this in every plotting script
# =============================================================================
#
# Based on Wong (2011, Nature Methods) + Tol colorblind-safe palettes.
# No red-green pairing. Perceptually distinct under deuteranopia & protanopia.
#
# Usage: source("colors_master.R") at the top of any plotting script.
#
# =============================================================================

# --- Diagnosis (2 levels) ---
# Wong blue vs Wong vermillion — maximally distinct under all CVD types
col_nd  <- "#D55E00"
col_t1d <- "#0072B2"
colors_diag <- c("ND" = col_nd, "T1D" = col_t1d)

# --- Region (3 levels) ---
# Wong amber / sky blue / bluish green — proven colorblind-safe triad
col_head <- "#E69F00"
col_body <- "#56B4E9"
col_tail <- "#009E73"
colors_region <- c("Head" = col_head, "Body" = col_body, "Tail" = col_tail)

# --- Sex (2 levels) ---
# Wong reddish purple / Tol dark indigo — distinct from Diag and Region
col_female <- "#CC79A7"
col_male   <- "#332288"
colors_sex <- c("Female" = col_female, "Male" = col_male)

# --- Disease duration bins (sequential, T1D only) ---
# Orange sequential scale — never appears alongside other categories
col_dur_short  <- "#FDAE61"  # 0-5 yr
col_dur_mid    <- "#F46D43"  # 5-15 yr
col_dur_long   <- "#A50026"  # >15 yr
colors_duration <- c("0\u20135 yr" = col_dur_short,
                     "5\u201315 yr" = col_dur_mid,
                     ">15 yr"      = col_dur_long)

# =============================================================================
# NOTES FOR REVIEWERS / COLLABORATORS
# =============================================================================
# 
# Category     Hex       Source          CVD-safe?
# ─────────    ───────   ─────────────   ─────────
# ND           #0072B2   Wong blue       Yes
# T1D          #D55E00   Wong vermillion Yes
# Head         #E69F00   Wong orange     Yes
# Body         #56B4E9   Wong sky blue   Yes
# Tail         #009E73   Wong bl. green  Yes
# Female       #CC79A7   Wong red-purple Yes
# Male         #332288   Tol indigo      Yes
# Dur 0-5      #FDAE61   Sequential      Yes (only w/ other durations)
# Dur 5-15     #F46D43   Sequential      Yes
# Dur >15      #A50026   Sequential      Yes
#
# Categories that share a color family (ND-blue ~ Body-skyblue;
# T1D-vermillion ~ Head-amber) never appear in the same legend.
# They are always separated by faceting.
# =============================================================================
