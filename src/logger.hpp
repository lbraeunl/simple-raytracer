#pragma once

inline constexpr bool HEATMAP = true;
inline int INT_COUNT = 0;

#define LOG(code) do { if constexpr (HEATMAP) { code; } } while(0)