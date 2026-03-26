# PRP_DimensionSwitch
Internal repo for Ashwin &amp; Atsushi

For storing scripts, figures, models, etc. 


Here is a detailed description of the **DimensionSwitch_PRP** cognitive psychology task: *(Thank you Claude Sonnet 4.6!)*

---

## Overview

This is a **Psychological Refractory Period (PRP)** paradigm combined with a **Dimension Switching** task. Participants must respond to two stimuli presented in rapid succession, each requiring a different classification judgment. The key manipulation is the **Stimulus Onset Asynchrony (SOA)** — the time gap between the two stimuli — to examine how response selection for one task interferes with response selection for another.

---

## Stimuli

Each trial presents **two bivariate stimuli** (two dimensions each), displayed simultaneously on screen — one above center and one below center. The two stimulus dimensions are:

**Dimension 1 — Number:** drawn from the set `{3, 4, 6, 7}`
**Dimension 2 — Color:** drawn from the set `{Red, Light Red, Light Blue, Blue}`

So a single stimulus might be, for example, a **"3" printed in blue**, or a **"7" printed in light red**. On every trial, two such stimuli appear: one positioned above the screen center and one below, with their positions randomized each trial.

The colors map to specific RGB values:
- **Red:** (255, 0, 0)
- **Light Red:** (250, 192, 203) — a soft pink
- **Blue:** (0, 0, 205) — medium blue
- **Light Blue:** (135, 206, 250) — sky blue

Stimuli are rendered in a large font (size 100) against a **gray background** (165, 165, 165), within white rectangular frames (200×200 px) offset ±150 pixels vertically from center.

---

## Task Rules

There are **four possible classification tasks**, organized into two families by stimulus dimension:

| Family | Task Name | Classification Rule |
|--------|-----------|---------------------|
| Number | **Low–High** | Is the number low (3,4) or high (6,7)? |
| Number | **Odd–Even** | Is the number odd (3,7) or even (4,6)? |
| Color | **Red–Blue** | Is the color red-ish or blue-ish? |
| Color | **Bold–Faint** | Is the color bold/saturated or faint/light? |

Each **block** uses a **pair of tasks** — one from each family (e.g., Low–High paired with Red–Blue). The task pairing is assigned per block and shown to participants at the start of each block in a block introduction screen (e.g., *"1: Low-High  2: Red-Blue"*). An instruction image is also displayed.

Because the two tasks always draw on different stimulus dimensions (one on number, one on color), Stimulus 1 is always categorized on one dimension and Stimulus 2 on the other. Crucially, each stimulus object contains *both* a number and a color, so the participant must attend selectively to the task-relevant dimension.

**Stimulus-Response Mapping:** Two response keys are used per task (left/right), assigned to the two categories. Keys are:
- Top stimulus response: **A** (left) / **S** (right)
- Bottom stimulus response: **Z** (left) / **X** (right)

The correct key for each stimulus is determined by the task rule and the stimulus-response reference table (`SRREF`), which maps each of the 4 stimulus values to one of two response categories.

---

## Trial Structure & Timing

Each trial proceeds through the following sequence:

### 1. RSI / Fixation (300 ms)
A **neutral white fixation cross (+)** is displayed at the center of the screen, flanked by two empty white rectangular frames (one above, one below center). The inter-trial interval (RSI) is **300 ms** (450 ms for the very first trial of a block).

### 2. Stimulus 1 Onset (time = 0)
The first stimulus appears inside one of the two frames (top or bottom, randomly assigned). It is displayed as a colored number.

### 3. Stimulus 2 Onset (time = SOA)
The second stimulus appears in the other frame. The **SOA** (Stimulus Onset Asynchrony) between the two onsets is randomly sampled from: **100, 250, 500, 750, or 1000 ms**.

### 4. Stimulus Offset
Each stimulus remains visible for **200 ms** from its onset, then disappears. So both stimuli eventually go off, but possibly at different times depending on the SOA.

### 5. Response Collection
The participant must press the correct key for **each** of the two stimuli. Two responses are required per trial; the trial loop continues until both responses are recorded (or the block time elapses). Responses are recorded with millisecond precision relative to each stimulus's onset.

### 6. Feedback (100 ms)
Immediately after both responses are collected, the fixation cross changes color:
- **Green** = both responses correct
- **Red** = one or both responses incorrect

The feedback is brief — **100 ms** — then the screen resets and the next trial's RSI begins.

---

## Block Structure

- **1 practice block** (~15 seconds) followed by **150 experimental blocks** (~20 seconds each)
- Within each block, trials repeat continuously until the block's time limit is reached, up to a maximum of **35 trials per block**
- The task pair is fixed within a block but varies across blocks (randomly assigned)
- At the end of each block, participants see a **results screen** showing how many trials they missed

---

## Incentives / Performance-Based Points

Participants can earn points each block based on three criteria:
1. **Accuracy** ≥ 80% correct
2. **At least 7 responses** completed in the block
3. **RT** below the 70th percentile of all prior experimental RTs (i.e., responding fast relative to their own history)

If all three criteria are met, points are awarded equal to the number of fast-RT trials (halved and rounded up). A running total of points is displayed at the end of each block.

---

## EEG / Neuroimaging Integration

The task includes full **EEG trigger coding** via a parallel port. Event codes are sent for:
- Block onset (`code 10`)
- Trial onset (`code 9`)
- Stimulus 1 onset (`code 1`) and Stimulus 2 onset (`code 2`)
- Response accuracy (`code 11` = correct, `code 10` = incorrect)
- Characteristic codes encoding block number, trial number, and full trial condition

This indicates the task was designed for use with **EEG** (and possibly ECoG) to examine neural correlates of dual-task response selection and task switching.

---

## Summary

In essence, this is a fast-paced dual-task paradigm where participants classify two overlapping colored-number stimuli, each by a different rule (e.g., is this number odd/even? is that color red/blue?). The SOA between the stimuli is the primary independent variable, allowing researchers to probe how the brain manages two competing response-selection processes — the classic signature of the **PRP bottleneck effect**. The added dimension-switching element (task pairs vary by block) allows investigation of how task-set reconfiguration interacts with dual-task costs.
