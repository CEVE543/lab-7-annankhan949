// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  pagenumbering: "1",
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: pagenumbering,
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(center)[#block(inset: 2em)[
      #set par(leading: heading-line-height)
      #if (heading-family != none or heading-weight != "bold" or heading-style != "normal"
           or heading-color != black or heading-decoration == "underline"
           or heading-background-color != none) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#import "@preview/fontawesome:0.5.0": *

#show: doc => article(
  title: [CEVE 543 Fall 2025 Lab 7: Bias Correction Implementation],
  subtitle: [Delta method and quantile mapping for temperature bias correction],
  authors: (
    ( name: [ak276],
      affiliation: [],
      email: [] ),
    ),
  date: [2025-10-24],
  margin: (x: 1in,y: 1in,),
  fontsize: 11pt,
  sectionnumbering: "1.1.a",
  pagenumbering: "1",
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)

= Background and Goals
<background-and-goals>
Climate models have systematic biases that directly affect impact assessments, arising from coarse spatial resolution, parameterization of sub-grid processes, representation of topography, and errors in simulated circulation patterns. This lab implements two widely-used bias correction approaches: the delta method and quantile-quantile (QQ) mapping. The delta method preserves the climate model's change signal while anchoring absolute values to observations, whereas QQ-mapping corrects the full distribution of values.

Both methods assume stationarity---that the statistical relationship between model and observations remains constant across climate states. This assumption may not hold under significant climate change. We'll explore the strengths and limitations of each method using temperature data for Boston, providing hands-on experience before PS2 Part 1.

= Study Location and Data
<study-location-and-data>
This lab uses temperature data for Boston Logan International Airport (Station ID: USW00014739, 42.3631$""^circle.stroked.tiny$N, 71.0064$""^circle.stroked.tiny$W). Observational data comes from GHCN-Daily (1936-2024) and is provided in the `USW00014739.csv` file in degrees Celsius. Climate model data comes from the GFDL-ESM4 model's 3-hourly near-surface air temperature (`tas`), pre-downloaded from Google Cloud Storage for both historical (1850-2014) and SSP3-7.0 (2015-2100) scenarios. Refer to #link("./labs/Lab-6/index.qmd")[Lab 6] for details on CMIP6 data structure.

== Data Processing Notes
<data-processing-notes>
The GHCN data provides daily average temperature in degrees Celsius. CMIP6 provides 3-hourly instantaneous temperature in Kelvin, which we'll convert to Celsius and aggregate to daily averages. The pre-downloaded NetCDF files (`boston_historical.nc` and `boston_ssp370.nc`) contain 3-hourly surface air temperature for the grid cell nearest to Boston. We aggregate 8 consecutive 3-hour periods to approximate daily averages, ignoring daylight saving time---a simplification typical in many bias correction applications. If you're interested in applying these methods to a different location, the `download_data.jl` script demonstrates how to extract CMIP6 data from Google Cloud Storage.

#block[
#callout(
body: 
[
Before starting the lab, uncomment the `Pkg.instantiate()` line in the first code block and run it to install all required packages. This will take a few minutes the first time. After installation completes, comment the line back out to avoid reinstalling on subsequent runs.

]
, 
title: 
[
Before Starting
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
#callout(
body: 
[
The render was not working because it was struggling to load bibliography("../../references.bib"). I got rid of it so the render would work, originally tried to grab reference.bib from older labs and add into directory but it still would not print.

]
, 
title: 
[
references.bib
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
= Lab Implementation
<lab-implementation>
== Package Setup
<package-setup>
#block[
```julia
using Pkg
lab_dir = dirname(@__FILE__)
Pkg.activate(lab_dir)
#Pkg.instantiate() # uncomment this the first time you run the lab to install packages, then comment it back
```

]
#block[
```julia
using CSV, CairoMakie, DataFrames, Dates, LaTeXStrings, NCDatasets, Statistics, StatsBase, TidierData
ENV["DATAFRAMES_ROWS"] = 6
CairoMakie.activate!()

# Constants for plotting
const MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
const MONTH_LABELS = (1:12, MONTH_NAMES)
```

]
== Task 1: Load and Process Data
<task-1-load-and-process-data>
We begin by loading both observational and climate model data. The observational data from GHCN-Daily spans 1936-2024, while the CMIP6 data requires processing from 3-hourly to daily resolution. We'll use the overlapping historical period (1995-2014) to calibrate bias corrections.

#block[
#callout(
body: 
[
Load the Boston GHCN temperature data from `USW00014739.csv` and filter to years with at least 80% complete data. Then load the pre-downloaded CMIP6 NetCDF files and aggregate the 3-hourly temperature data to daily values for the historical (1995-2014) and near-term future (2020-2040) periods. The helper functions `load_cmip6_data()` and `aggregate_to_daily()` are provided below. Finally, visualize the annual cycle comparing observations vs the historical GCM simulation.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
# Load and clean the observational data
data_path = joinpath(lab_dir, "USW00014739.csv")
df = @chain begin
    CSV.read(data_path, DataFrame)
    @mutate(
        TAVG = ifelse.(ismissing.(TMIN) .| ismissing.(TMAX), missing, (TMIN .+ TMAX) ./ 2),
    )
    @mutate(TAVG = TAVG / 10.0) # Convert to degrees C
    @rename(date = DATE)
    @mutate(year = year(date), month = month(date))
    @select(date, year, month, TAVG)
end

# Filter to years with at least 80% complete data
yearly_counts = @chain df begin
    @group_by(year)
    @summarize(n_obs = sum(!ismissing(TAVG)), n_total = n())
    @mutate(frac_complete = n_obs / n_total)
end

good_years_df = @chain yearly_counts begin
    @filter(frac_complete >= 0.8)
end
good_years = good_years_df.year

df_clean = @chain df begin
    @filter(year in !!good_years)
    dropmissing(:TAVG)
end
```

]
#block[
```julia
# Helper functions for CMIP6 data
"""
Load CMIP6 temperature data from a local NetCDF file and convert times to DateTime.
"""
function load_cmip6_data(file_path::String)
    ds = NCDataset(file_path)
    tas_data = ds["tas"][:]
    time_cf = ds["time"][:]
    close(ds)

    time_data = [DateTime(
        Dates.year(t), Dates.month(t), Dates.day(t),
        Dates.hour(t), Dates.minute(t), Dates.second(t)
    ) for t in time_cf]

    return tas_data, time_data
end

"""
Aggregate 3-hourly temperature data to daily averages.
"""
function aggregate_to_daily(tas_3hr, time_3hr)
    n_3hr_per_day = 8
    n_days = div(length(tas_3hr), n_3hr_per_day)

    daily_temp = Vector{Float64}(undef, n_days)
    daily_dates = Vector{Date}(undef, n_days)

    for i in 1:n_days
        idx_start = (i - 1) * n_3hr_per_day + 1
        idx_end = i * n_3hr_per_day
        daily_vals = collect(skipmissing(tas_3hr[idx_start:idx_end]))
        daily_temp[i] = isempty(daily_vals) ? NaN : mean(daily_vals)
        daily_dates[i] = Date(time_3hr[idx_start])
    end

    daily_temp_c = daily_temp .- 273.15  # Convert K to C
    return daily_temp_c, daily_dates
end

# Load and process CMIP6 data
hist_file = joinpath(lab_dir, "boston_historical.nc")
ssp370_file = joinpath(lab_dir, "boston_ssp370.nc")

tas_hist_3hr, time_hist_3hr = load_cmip6_data(hist_file)
tas_ssp370_3hr, time_ssp370_3hr = load_cmip6_data(ssp370_file)

# Process historical period (1995-2014)
hist_start = DateTime(1995, 1, 1)
hist_end = DateTime(2014, 12, 31, 23, 59, 59)
hist_idx = (time_hist_3hr .>= hist_start) .& (time_hist_3hr .<= hist_end)
tas_hist_daily, dates_hist_daily = aggregate_to_daily(tas_hist_3hr[hist_idx], time_hist_3hr[hist_idx])

# Process near-term period (2020-2040)
near_start = DateTime(2020, 1, 1)
near_end = DateTime(2040, 12, 31, 23, 59, 59)
near_idx = (time_ssp370_3hr .>= near_start) .& (time_ssp370_3hr .<= near_end)
tas_ssp370_near_daily, dates_ssp370_near_daily = aggregate_to_daily(tas_ssp370_3hr[near_idx], time_ssp370_3hr[near_idx])

# Create DataFrames
df_gcm_hist = DataFrame(
    date=dates_hist_daily, temp=tas_hist_daily,
    year=year.(dates_hist_daily), month=month.(dates_hist_daily)
)

df_ssp370_near = DataFrame(
    date=dates_ssp370_near_daily, temp=tas_ssp370_near_daily,
    year=year.(dates_ssp370_near_daily), month=month.(dates_ssp370_near_daily)
)

df_obs_hist = @chain df_clean begin
    @filter(year >= 1995 && year <= 2014)
end

# Compute monthly climatologies
obs_monthly = @chain df_obs_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(TAVG))
end

gcm_monthly = @chain df_gcm_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(temp))
end
```

]
```julia
let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel="Month",
        ylabel=L"Temperature ($^\circ$C)",
        title="Annual Cycle: Observations vs GCM Historical (1995-2014)",
        xticks=MONTH_LABELS)
    lines!(ax, obs_monthly.month, obs_monthly.mean_temp, linewidth=2, color=:steelblue, label="Observations")
    lines!(ax, gcm_monthly.month, gcm_monthly.mean_temp, linewidth=2, color=:coral, label="GCM Historical")
    axislegend(ax, position=:lt)
    fig
end
```

#figure([
#box(image("index_files/figure-typst/fig-gcm-obs-comparison-output-1.png", height: 3.5in, width: 5.5in))
], caption: figure.caption(
position: bottom, 
[
Annual cycle comparison between GHCN observations and GFDL-ESM4 historical simulation for Boston (1995-2014). The GCM shows a warm bias across most months.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-gcm-obs-comparison>


== Task 2: Implement Delta Method
<task-2-implement-delta-method>
The delta method corrects the mean bias while preserving the climate model's projected change signal. We calculate a monthly bias correction based on the historical period (1995-2014) and apply it to future projections.

#block[
#callout(
body: 
[
Implement the additive delta method for temperature bias correction. For each calendar month $m$, calculate the mean bias: $Delta_m = macron(T)_(upright("hist") \, m)^(upright("GCM")) - macron(T)_(upright("hist") \, m)^(upright("obs"))$. Then apply the correction to future values: $T_(upright("fut"))^(upright("corr")) (d \, m \, y) = T_(upright("fut"))^(upright("GCM")) (d \, m \, y) - Delta_m$.

Follow these steps:

+ Calculate the monthly mean bias by grouping both `df_gcm_hist` and `df_obs_hist` by month, computing their means, joining them, and computing the difference.
+ Create a bar plot visualizing the monthly bias.
+ Write a function `apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)` that applies the bias correction to a vector of temperatures and dates.
+ Apply your function to the near-term data and add a new column `temp_delta` to `df_ssp370_near`.
+ Create monthly climatologies for both raw and delta-corrected near-term data.
+ Visualize the annual cycle showing historical observations, raw GCM near-term, and delta-corrected near-term temperatures.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia

# 1. Compute monthly mean bias
monthly_bias = @chain begin
    df_gcm_hist
    @group_by(month)
    @summarize(gcm_mean = mean(temp))
    @left_join(obs_monthly, month)
    @mutate(bias = gcm_mean - mean_temp)
end
```

```
12×4 DataFrame
 Row │ month  gcm_mean  mean_temp   bias     
     │ Int64  Float64   Float64?    Float64  
─────┼───────────────────────────────────────
   1 │     1  -3.18274  -1.09024    -2.0925
   2 │     2  -1.70476   0.0218584  -1.72662
   3 │     3   2.00539   3.78944    -1.78405
   4 │     4   6.41217   9.34408    -2.93191
   5 │     5  12.445    14.5874     -2.14243
   6 │     6  17.8571   19.9404     -2.08332
   7 │     7  21.1401   23.3412     -2.20107
   8 │     8  20.5512   22.5353     -1.9841
   9 │     9  17.1836   18.6946     -1.51102
  10 │    10  11.555    12.6931     -1.1381
  11 │    11   5.14819   7.06367    -1.91548
  12 │    12   0.52005   2.04863    -1.52858
```

```julia
# 2. Visualize the bias
let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel = "Month",
        ylabel = "Temperature Bias (Celsius)",
        title = "Monthly Mean Temperature Bias (GCM - Observations)",
        xticks = MONTH_LABELS)
    barplot!(ax, monthly_bias.month, monthly_bias.bias, color = :teal)
    fig
end
```

#box(image("index_files/figure-typst/cell-8-output-1.png", height: 3.5in, width: 5.5in))

```julia
# 3. Defining apply_delta_method
function apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)
    temps_withBias = similar(gcm_temps)
    for i in eachindex(gcm_temps)
        m = month(gcm_dates[i]) 
        bias_val = monthly_bias_df[monthly_bias_df.month .== m, :bias][1]# Look up monthly bias
        temps_withBias[i] = gcm_temps[i] - bias_val
    end
    return temps_withBias
end
```

```
apply_delta_method (generic function with 1 method)
```

```julia
# 4. Apply delta correction to near-term GCM data, add in Delta Temp to df
df_ssp370_near.DeltaTemp = apply_delta_method(df_ssp370_near.temp, df_ssp370_near.date, monthly_bias)
```

```
7665-element Vector{Float64}:
 -5.293100300450412
 -5.623239460606662
 -0.9786776930285366
  2.5020046800183384
  5.411062297205838
  0.4519863694714634
 -3.4730624586535366
 -6.282297077794162
  1.1962185472058384
 -2.8390292555285366
  ⋮
  2.642983674080141
  0.09229397681451634
 -0.6736667165448587
  0.44959378150201634
 -1.1007602224042337
 -1.8960483083417337
  1.2649624338457663
 -1.7724215993573587
 -3.997549773185484
```

```julia
# 5. Compute monthly climatologies for raw and corrected GCM data
near_monthly_raw = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(temp))
end

near_monthly_delta = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(DeltaTemp))
end
```

```
12×2 DataFrame
 Row │ month  mean_temp 
     │ Int64  Float64   
─────┼──────────────────
   1 │     1   0.103505
   2 │     2   1.70296
   3 │     3   4.62097
   4 │     4  10.7028
   5 │     5  14.9487
   6 │     6  20.0121
   7 │     7  23.8329
   8 │     8  23.6061
   9 │     9  19.4281
  10 │    10  13.4868
  11 │    11   7.91067
  12 │    12   2.51416
```

```julia
# 6. Comparison plot: Observations vs GCM Raw vs Delta-Corrected
let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel = "Month",
        ylabel = "Temperature (Celsius)",
        title = "Annual Cycle for SSP3-7.0 Near-Term (2020–2040)",
        xticks = MONTH_LABELS)
    lines!(ax, obs_monthly.month, obs_monthly.mean_temp, color = :black, linewidth = 2, label = "Historical Obs")
    lines!(ax, near_monthly_raw.month, near_monthly_raw.mean_temp, color = :red, linewidth = 2, label = "GCM Raw")
    lines!(ax, near_monthly_delta.month, near_monthly_delta.mean_temp,
        color = :teal, linewidth = 2, label = "Delta Corrected")
    axislegend(ax, position = :lt)
    fig
end
```

#box(image("index_files/figure-typst/cell-12-output-1.png", height: 3.5in, width: 5.5in))

== Task 3: Implement Quantile-Quantile Mapping
<task-3-implement-quantile-quantile-mapping>
Unlike the delta method, QQ-mapping transforms the entire probability distribution of model output to match observations. For each value in the future model output, we find its percentile in the historical model distribution, then map it to the same percentile in the historical observed distribution.

#block[
#callout(
body: 
[
Implement QQ-mapping to correct the full distribution of temperature values.

Follow these steps:

+ Extract the observed and GCM historical temperatures as vectors (`obs_hist_temps` and `gcm_hist_temps`).
+ Use `ecdf()` from StatsBase.jl to fit empirical cumulative distribution functions to both datasets.
+ Write a function `apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)` that:
  - For each temperature value, finds its percentile in the GCM CDF
  - Clamps the percentile to the range \[0.001, 0.999\] to avoid extrapolation
  - Maps to the same percentile in the observed distribution using `quantile()`
+ Apply your function to the near-term data and add a new column `temp_qqmap` to `df_ssp370_near`.
+ Create a QQ-plot comparing observed and GCM quantiles for the historical period, showing the 1:1 line.

#block[
#callout(
body: 
[
- Use `ecdf_gcm(value)` to get the percentile of a value in the GCM distribution
- Use `quantile(obs_hist_temps, p)` to get the value at percentile `p` in the observed distribution
- The `clamp(x, low, high)` function constrains `x` to the range \[low, high\]

]
, 
title: 
[
Hints
]
, 
background_color: 
rgb("#ccf1e3")
, 
icon_color: 
rgb("#00A047")
, 
icon: 
fa-lightbulb()
, 
body_background_color: 
white
)
]
]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
using StatsBase
```

]
```julia
# 1. Extract historical temperatures as vectors
show(df_obs_hist) #TAVG
show(df_gcm_hist) #temp

obs_hist_temps = collect(df_obs_hist.TAVG)
gcm_hist_temps = collect(df_gcm_hist.temp)
```

#block[
```
7305×4 DataFrame
  Row │ date        year   month  TAVG    
      │ Date        Int64  Int64  Float64 
──────┼───────────────────────────────────
    1 │ 1995-01-01   1995      1     3.9
    2 │ 1995-01-02   1995      1     1.35
    3 │ 1995-01-03   1995      1    -0.8
    4 │ 1995-01-04   1995      1    -3.05
    5 │ 1995-01-05   1995      1    -6.35
    6 │ 1995-01-06   1995      1    -0.8
    7 │ 1995-01-07   1995      1     6.95
    8 │ 1995-01-08   1995      1    -0.55
  ⋮   │     ⋮         ⋮      ⋮       ⋮
 7299 │ 2014-12-25   2014     12    10.85
 7300 │ 2014-12-26   2014     12     6.4
 7301 │ 2014-12-27   2014     12     7.2
 7302 │ 2014-12-28   2014     12     7.75
 7303 │ 2014-12-29   2014     12     2.0
 7304 │ 2014-12-30   2014     12    -3.25
 7305 │ 2014-12-31   2014     12    -3.8
                         7290 rows omitted7300×4 DataFrame
  Row │ date        temp        year   month 
      │ Date        Float64     Int64  Int64 
──────┼──────────────────────────────────────
    1 │ 1995-01-01   -2.85105    1995      1
    2 │ 1995-01-02    1.71316    1995      1
    3 │ 1995-01-03    7.35714    1995      1
    4 │ 1995-01-04   -2.98737    1995      1
    5 │ 1995-01-05   -3.1699     1995      1
    6 │ 1995-01-06   -9.43973    1995      1
    7 │ 1995-01-07   -8.61146    1995      1
    8 │ 1995-01-08   -7.47367    1995      1
  ⋮   │     ⋮           ⋮         ⋮      ⋮
 7294 │ 2014-12-25    0.607141   2014     12
 7295 │ 2014-12-26   -4.0741     2014     12
 7296 │ 2014-12-27   -3.99879    2014     12
 7297 │ 2014-12-28    5.19061    2014     12
 7298 │ 2014-12-29    1.59582    2014     12
 7299 │ 2014-12-30   -0.54917    2014     12
 7300 │ 2014-12-31    7.50701    2014     12
                            7285 rows omitted
```

]
```
7300-element Vector{Float64}:
  -2.8510498046874773
   1.7131591796875227
   7.357141113281273
  -2.9873718261718523
  -3.1698974609374773
  -9.439733886718727
  -8.611456298828102
  -7.473669433593727
 -13.736547851562477
 -10.071051025390602
   ⋮
  -1.3686889648437273
   9.637475585937523
   0.6071411132812727
  -4.074102783203102
  -3.9987854003906023
   5.190606689453148
   1.5958190917968977
  -0.5491699218749773
   7.507012939453148
```

```julia
# 2. Fit empirical CDFs
ecdf_gcm = ecdf(gcm_hist_temps)
ecdf_obs = ecdf(obs_hist_temps)
```

```
ECDF{Vector{Float64}, Weights{Float64, Float64, Vector{Float64}}}([-16.95, -16.65, -15.85, -15.0, -14.75, -14.7, -13.9, -13.6, -13.35, -13.35  …  30.6, 30.85, 30.85, 31.1, 31.15, 31.15, 31.35, 31.4, 31.95, 33.3], Float64[])
```

```julia
# 3. Define QQ-mapping function
function apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)
    QQ_temps = similar(gcm_temps)
    for i in eachindex(gcm_temps)
        p = ecdf_gcm(gcm_temps[i])# Find percentile using `ecdf_gcm(value)` 
        p_clamp = clamp(p, 0.001, 0.999)#The `clamp(x, low, high)` function constrains `x` to the range [low, high]
        QQ_temps[i] = quantile(obs_hist_temps, p_clamp)# Map to observed quantile. Using `quantile(obs_hist_temps, p)`
    end
    return QQ_temps
end
```

```
apply_qqmap_empirical (generic function with 1 method)
```

```julia
# 4. Apply QQ-mapping correction to near-term (2020–2040) GCM data
df_ssp370_near.QQTemp = apply_qqmap_empirical(df_ssp370_near.temp, ecdf_gcm, obs_hist_temps)
```

```
7665-element Vector{Float64}:
 -5.75
 -6.1
 -1.95
  2.5
  5.8
 -0.3
 -4.15
 -6.65
  0.8
 -3.65
  ⋮
  3.35
  0.0
 -1.1
  0.55
 -1.65
 -2.25
  1.65
 -2.2
 -4.15
```

```julia
# 5. Compute quantiles at evenly spaced probabilities and plot 

probs = 0.01:0.01:0.99 #Range of percentiles from 1% to 99% in 1% increments
obs_quantiles = [quantile(obs_hist_temps, p) for p in probs]
gcm_quantiles = [quantile(gcm_hist_temps, p) for p in probs]

let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel = "GCM Historical Quantiles (Celsius)",
        ylabel = "Observed Historical Quantiles (Celsius)",
        title = "GCM vs Observed Historical (QQ Plot 1995–2014)")
    scatter!(ax, gcm_quantiles, obs_quantiles, color = :black, markersize = 6, label = "GCM Quantiles")
    lines!(ax, gcm_quantiles, gcm_quantiles, color = :red, label = "1:1 Line")
    axislegend(ax, position = :lt)
    fig
end
```

#box(image("index_files/figure-typst/cell-18-output-1.png", height: 3.5in, width: 5.5in))

== Task 4: Compare Methods
<task-4-compare-methods>
Now we compare the delta method and QQ-mapping approaches to understand their strengths and limitations.

#block[
#callout(
body: 
[
Create a single figure showing monthly mean temperature for the near-term period (2020-2040) using all four approaches:

+ Compute the monthly mean for QQ-mapped temperatures (`near_monthly_qqmap`)
+ Create a figure with 4 lines:
  - Historical observations (`obs_monthly`)
  - Raw GCM near-term (`near_monthly_raw`)
  - Delta-corrected near-term (`near_monthly_delta`)
  - QQ-mapped near-term (`near_monthly_qqmap`)

After creating the plot, examine it carefully and consider:

- How do the two correction methods differ?
- What does the QQ-mapped line reveal about this method's treatment of the warming signal?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia
#Monthly mean for QQ
near_monthly_qqmap = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(QQTemp))
end

#Plot the 4 lines
let
    fig = Figure(resolution = (800, 450))
    ax = Axis(fig[1, 1],
        xlabel = "Month",
        ylabel = "Temperature (Celsius)",
        title = "Near-Term (2020–2040) Monthly Means: Obs vs GCM (Raw, Delta, QQ-mapped)",
        xticks = MONTH_LABELS)
    lines!(ax, obs_monthly.month, obs_monthly.mean_temp,linewidth = 2, label = "Historical Obs", linestyle = :solid)
    lines!(ax, near_monthly_raw.month, near_monthly_raw.mean_temp, linewidth = 2, label = "GCM Raw", linestyle = :solid)
    lines!(ax, near_monthly_delta.month, near_monthly_delta.mean_temp,
        linewidth = 2, label = "Delta Corrected", linestyle = :dash)
    lines!(ax, near_monthly_qqmap.month, near_monthly_qqmap.mean_temp,
        linewidth = 2, label = "QQ-mapped", linestyle = :dot)
    axislegend(ax, position = :lt)
    fig
end
```

#block[
```
┌ Warning: Found `resolution` in the theme when creating a `Scene`. The `resolution` keyword for `Scene`s and `Figure`s has been deprecated. Use `Figure(; size = ...` or `Scene(; size = ...)` instead, which better reflects that this is a unitless size and not a pixel resolution. The key could also come from `set_theme!` calls or related theming functions.
└ @ Makie ~/.julia/packages/Makie/4JW9B/src/scenes.jl:264
```

]
#box(image("index_files/figure-typst/cell-19-output-2.png", height: 4.6875in, width: 8.33333in))

== Task 5: Reflection
<task-5-reflection>
Finally, we reflect on the assumptions, limitations, and appropriate use cases for each bias correction method.

#block[
#callout(
body: 
[
Write brief responses (2-3 sentences each) to the following questions:

+ #strong[Method selection:] If you were providing climate data to support a decision about urban heat management in Boston, which bias correction method would you recommend and why?

+ #strong[Appropriateness of QQ-mapping:] In the #cite(<ines_biascorrection:2006>, form: "prose") paper we discussed, QQ-mapping was used to correct both the #emph[frequency] and #emph[intensity] of rainfall for crop modeling. For temperature data, does it make sense to correct the full distribution? Why or why not? When might the delta method be more appropriate than QQ-mapping?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
=== Method Selection for Urban Heat Management
<method-selection-for-urban-heat-management>
#emph[I would reccomend QQ Mapping as the bias correction method. With using a delta method it assumes that the bias is constant over time, which makes it less viable at factoring in extremes like heat. QQ Mapping allows it to be able to adjusted based on percentiles rather than just shifting the mean like with the delta method.]

=== Appropriateness of QQ-Mapping for Temperature
<appropriateness-of-qq-mapping-for-temperature>
#emph[In this case study, yes I believe it is appropriate to use QQ Mapping. Unless we are looking at large scale temperature shifts, I do not believe that the bias method would be as effective. For studies regarding heat management, where there is a need to identify extremes, QQ mapping is more appropriate.]

= References
<references>
#block[
] <refs>




