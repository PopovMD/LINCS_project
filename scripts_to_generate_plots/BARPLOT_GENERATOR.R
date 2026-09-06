#!/usr/bin/env rscript

################### BARPLOT_GENERATOR.R #########################
#
# args[1] --- name of inhibitor's target (MEK, mTOR, PI3K, ...)
# args[2] --- name of an algorithm (L2S2, iLINCS, L1000, ...)
# args[3] --- style of plots ("dark" or "white")
# args[4] --- pathway from root to the projects' folder
#
#################################################################


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(jsonlite)



args <- commandArgs(trailingOnly = TRUE)

inh_grp = args[1] #"PI3K"
method = args[2] #"iLINCS"
plots_mode = args[3] #"wht"
project_folder_name = args[4] #"/Users/svetlana/Desktop/CMap_2025-26/folder_to_git/"

if (method == "L2S2") {
  res_dir_suff = "_results_20000"
} else if (method == "iLINCS") {
  res_dir_suff = "_results"
} else if (method == "L1000") {
  res_dir_suff = "_results_with_empty_files"
}
if (inh_grp == "MEK") {
  annotation_segment <- c(1:41)
} else if (inh_grp == "mTOR") {
  annotation_segment <- c(42:135)
} else if (inh_grp == "PI3K") {
  annotation_segment <- c(140:148, 150:154)
} 

method_res_dir = paste0(project_folder_name, inh_grp, "_inhibitors/", method, "/", method, res_dir_suff)
methodres <- list.files(method_res_dir)
methodres <- methodres[order(as.numeric(str_extract(methodres, "^\\d+")))]

annotation <- read_tsv(paste0(project_folder_name, "annotations_parsed/sgn-perturbagen_annotation.tsv"), col_names = T)
nrow(annotation)
synonyms <- fromJSON(paste0(project_folder_name, inh_grp, "_inhibitors/", inh_grp, "_perturbagens_synonyms.json"))


num_of_row_for_each_sgn <- vector(mode="numeric", length = length(methodres))
first_perts <- vector(mode="character", length = length(methodres))
ann_row_nums <- vector(mode="numeric", length = length(methodres))

nrow(annotation)



# заранее определяем названия колонок, по которым будем сортировать таблицы с результатами работы алгоритмов
if (method == "L2S2") {
  pval_col_name <- "adjPvalue"
  score_col_name <- "oddsRatio"
  pert_col_name <- "perturbation"
} else if (method == "iLINCS") {
  pval_col_name <- "pValue"
  score_col_name <- "similarity"
  pert_col_name <- "compound"
} else if (method == "L1000") {
  pval_col_name <- "p-value"
  score_col_name <- "Scores"
  pert_col_name <- "Perturbation"
}

setwd(method_res_dir)
###### СОБИРАЕМ НОМЕРА СТРОК ПЕРВЫХ ВХОЖДЕНИЙ ИСКОМОГО ВЕЩЕСТВА ######
for (i in 1:length(methodres)) {
  
  # считаем номер строки в annotation, которая соответствует текущей сигнатуре
  row_num = as.numeric(sub("^(\\d+).*", "\\1", methodres[i])) 
  ann_row_nums[i] <- row_num
  
  pert <- annotation[[row_num,2]]
  print(paste("PROCESSING SYMBOL:", methodres[i]))
  methodsgni <- read_tsv(methodres[i], show_col_types = FALSE)
  if (nrow(methodsgni) == 0) {
    num_of_row_for_each_sgn[i] <- NA
    next
  }
  
  ###### СОРТИРОВКА ######
  methodsgni[["Combined_score"]] <- (1 - methodsgni[[pval_col_name]]) * methodsgni[[score_col_name]]
  methodsgni <- methodsgni %>% arrange(desc(Combined_score))
  
  methodsgni <- methodsgni[[pert_col_name]]
  try(
    syns <- synonyms[[pert]]
  )
  
  # ищем первое вхождение какого-либо синонима (считаем ДОЛЮ)
  n <- which(toupper(methodsgni) %in% toupper(syns))
  if (length(n) == 0) {
    num_of_row_for_each_sgn[i] <- NA
    next
  }
  num_of_row_for_each_sgn[i] <- n[1] / length(methodsgni)
  try(first_perts[i] <- methodsgni[[n[1]]])
  
}

### СРАВНЕНИЕ РАЗНЫЙ АЙДИШНИКОВ И ПРОВЕРКА, ЧТО НИЧЕГО НЕ ПОТЕРЯЛИ
length(num_of_row_for_each_sgn)

num_of_row_for_each_sgn
sum(!(is.na(num_of_row_for_each_sgn)))
sum(!(is.na(num_of_row_for_each_sgn))) / length(num_of_row_for_each_sgn)
first_perts
ann_row_nums

annotation <- annotation %>% dplyr::slice(ann_row_nums)
annotation[["method_res"]] <- num_of_row_for_each_sgn
#annotation[["fp"]] <- first_perts
annotation$label_symbol <- paste0(round(annotation$method_res, 3), ", ", substr(annotation$perturbation, 1, 1))



## собираем ДЛИНЫ всех сигнатур ##
setwd(paste0(project_folder_name, inh_grp, "_inhibitors/signatures"))
sgns <- list.files()
sgns <- sgns[order(as.numeric(str_extract(sgns, "^\\d+")))]
sgn_lengths <- vector(mode="numeric", length = length(methodres))
#length(ilincsres) == length(sgns)

for (i in 1:length(sgns)) {
  print(paste("processing file №", i))
  f <- read_tsv(sgns[i], show_col_types = FALSE)
  sgn_lengths[i] <- nrow(f)
}

# теперь собираем индексы сигнатур, для которых есть результаты l1000
if (inh_grp == "PI3K") {
  sgn_coord_num <- 1:length(sgns)
} else {
  sgn_coord_num <- ann_row_nums - annotation_segment[1] + 1
}
annotation["sgn_lengths"] <- sgn_lengths[sgn_coord_num]



##### СОБИРАЕМ КОЛИЧЕСТВА ГЕНОВ В ПЕРЕCЕЧЕНИИ С LINCS в переменную PORTIONS_PER_SGN #####
lincs <- read_tsv(paste0(project_folder_name, "GSE92742_Broad_LINCS_gene_info.txt"))

tryCatch({
  setwd(paste0(project_folder_name, inh_grp, "_inhibitors/signatures"))
}, error = function(e) {
  print(paste0("No signatures` files for ", inh_grp, " inhibitors"))
})

sgns <- list.files(pattern = "*.txt")
# обязательно СОРТИРУЕМ сигнатуры по порядковым номерам
sgns <- sgns[order(as.numeric(str_extract(sgns, "^\\d+")))]

lm <- lincs %>%
  filter(pr_is_lm == 1) %>%
  pull(pr_gene_symbol) %>%
  unique()
print(paste("N of LM:", length(lm)))

all_lincs <- lincs %>%
  pull(pr_gene_symbol) %>%
  unique()
print(paste("N of all lincs genes:", length(all_lincs)))

all_genes <- c()
portions_per_sgn <- list(
  sgn_lm = c(),
  sgn_bing = c(),
  sgn_notlincs = c()
)


length(sgns)

for (file in sgns) {
  print(paste("PROCESSING", file))
  
  sgn <- tryCatch({
    read_tsv(file, show_col_types = FALSE) %>%
      pull(1) %>%
      unique()
  }, error = function(e) {
    print(paste("file", file, "is empty"))
    return(NULL)
  })
  
  if (is.null(sgn)) {
    portions_per_sgn$sgn_lm <- c(portions_per_sgn$sgn_lm, 0)
    portions_per_sgn$sgn_bing <- c(portions_per_sgn$sgn_bing, 0)
    portions_per_sgn$sgn_notlincs <- c(portions_per_sgn$sgn_notlincs, 0)
    next
  }
  
  all_genes <- unique(c(all_genes, sgn))
  
  sgn_lm <- length(intersect(sgn, lm))
  sgn_bing <- length(intersect(sgn, all_lincs)) - sgn_lm
  sgn_notlincs <- length(sgn) - sgn_bing - sgn_lm
  
  portions_per_sgn$sgn_lm <- c(portions_per_sgn$sgn_lm, sgn_lm)
  portions_per_sgn$sgn_bing <- c(portions_per_sgn$sgn_bing, sgn_bing)
  portions_per_sgn$sgn_notlincs <- c(portions_per_sgn$sgn_notlincs, sgn_notlincs)
}

length(portions_per_sgn$sgn_notlincs)

df <- tibble("sgn_num" = annotation_segment,
             "sgn_lm" = portions_per_sgn$sgn_lm,
             "sgn_bing" = portions_per_sgn$sgn_bing,
             "sgn_notlincs" = portions_per_sgn$sgn_notlincs)
df <- df %>% 
  pivot_longer(cols = c(sgn_lm, sgn_bing, sgn_notlincs),
               names_to = "gene_group",
               values_to = "gene_amounts")

distr <- ggplot(df, aes(x = sgn_num, y = gene_amounts, fill = gene_group)) +
  # "fill" растягивает столбцы до 100% (от 0 до 1)
  geom_col(position = "fill") + 
  scale_x_continuous(breaks = df$sgn_num) +
  theme(
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               vjust = 0.5,
                               size = 6)) +
  labs(y = "Доля", x = "Группа")

distr


all(portions_per_sgn$sgn_bing >= portions_per_sgn$sgn_lm)

annotation["n_of_lm"] <- portions_per_sgn$sgn_lm[sgn_coord_num]#[ann_row_nums]
annotation["n_of_bing"] <- portions_per_sgn$sgn_bing[sgn_coord_num]#[ann_row_nums]
annotation["n_of_notlincs"] <- portions_per_sgn$sgn_notlincs[sgn_coord_num]#[ann_row_nums]


##### CHECKING CORRELATIONS #####
corrdf <- annotation[!(is.na(annotation$method_res)), ]

res_len <- cor(corrdf$method_res, corrdf$sgn_lengths, method = "spearman")
res_lm <- cor(corrdf$method_res, corrdf$n_of_lm / corrdf$sgn_lengths, method = "spearman")
res_bing <- cor(corrdf$method_res, corrdf$n_of_bing/ corrdf$sgn_lengths, method = "spearman")
res_nl <- cor(corrdf$method_res, corrdf$n_of_notlincs / corrdf$sgn_lengths, method = "spearman")
res_lincs <- cor(corrdf$method_res, 
                 (corrdf$n_of_lm + corrdf$n_of_bing) / corrdf$sgn_lengths, 
                 method = "spearman")
res_lincs_abs <- cor(corrdf$method_res, 
                 (corrdf$n_of_lm + corrdf$n_of_bing), 
                 method = "spearman")
res_lm_abs <- cor(corrdf$method_res, corrdf$n_of_lm, method = "spearman")

table_new <- data.frame("res_len" = res_len, 
                        "res_lm" = res_lm, 
                        "res_bing" = res_bing,
                        "res_notLINCS" = res_nl,
                        "res_LINCS" = res_lincs,
                        "res_LINCS_abs" = res_lincs_abs)
table_new


##### BARPLOT FOR SYMBOLS#####
data_sorted_pert <- annotation %>%
  mutate(
    group = as.numeric(factor(perturbation)),
  ) %>%
  # Сортируем данные по этому новому числовому столбцу (чтоб сгруппировать по веществам)
  arrange(group)

# === НАЧАЛО ДОБАВЛЕННОГО КОДА ===
# Пересортировываем внутри каждой группы веществ по возрастанию длины сигнатуры
data_sorted_pert <- data_sorted_pert %>%
  arrange(group, sgn_lengths) %>% 
  filter(sgn_lengths != 0)
# Создаем фактор с уровнями в порядке сортировки данных
data_sorted_pert$sgn_number_ordered <- factor(data_sorted_pert$sgn_number, 
                                              levels = unique(data_sorted_pert$sgn_number))
# === КОНЕЦ ДОБАВЛЕННОГО КОДА ===

# ОТДЕЛЬНЫЙ БАРПЛОТ ДЛЯ НУМЕРОВ
bp_strnum <- ggplot(data = data_sorted_pert,
                    mapping = aes(sgn_number_ordered, 
                                  as.numeric(method_res),
                                  fill = perturbation)) +
  geom_col() +
  theme_minimal() +
  geom_text(
    aes(label = label_symbol), # Указываем, какой столбец использовать как текст метки
    hjust = 0,            # Вертикальное выравнивание: -0.5 немного приподнимает текст над столбцом
    color = "black",         # Цвет текста
    size = 2.5,                 # Размер текста
    angle = 90
  ) +
  theme(
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               vjust = 0.5,
                               size = 6)) +
  #scale_x_discrete(labels = function(x) str_replace(x, ".L2S2.tsv", ""))+
  labs(x = "Signatures", 
       y = paste0("Number of perturbagen in ", method, " result (Symbol)"), 
       title = paste0("Number of perturbagen in ", method, "\nresults per signature"))


bp_strnum



### преобразуем, чтобы нули корректно отображались на логарифмической кривой (нули остаются нулями, единицы - двойками и тд)
data_sorted_pert$sgn_lengths_transformed <- -log10(as.numeric(data_sorted_pert$sgn_lengths)+1)  # отрицательные значения для инверсии
data_sorted_pert$n_of_bing_transformed <- -log10(as.numeric(data_sorted_pert$n_of_bing)+1)
data_sorted_pert$n_of_lm_transformed <- -log10(as.numeric(data_sorted_pert$n_of_lm)+1)


# Подготовка данных
data_sorted_pert <- data_sorted_pert %>%
  mutate(
    # Исходные значения
    sgn_lengths_orig = as.numeric(sgn_lengths),
    n_of_bing_orig = as.numeric(n_of_bing),
    n_of_lm_orig = as.numeric(n_of_lm),
    n_of_notlincs_orig = as.numeric(n_of_notlincs),
    
    # Логарифмическая высота всего столбца
    total_height_log = -log10(sgn_lengths_orig + 1),
    
    # Вычисляем доли (на основе исходных, не логарифмированных значений)
    bing_prop = n_of_bing_orig / sgn_lengths_orig,
    lm_prop = n_of_lm_orig / sgn_lengths_orig,
    notlincs_prop = n_of_notlincs_orig / sgn_lengths_orig,
    
    # Вычисляем высоту каждого сегмента в логарифмическом пространстве
    bing_height_log = total_height_log * bing_prop,
    lm_height_log = total_height_log * lm_prop,
    notlincs_height_log = total_height_log * notlincs_prop,
    
    # Создаем фактор для правильного порядка
    sgn_number_factor = sgn_number_ordered #reorder(sgn_number, group)
  )

# Преобразуем в длинный формат для ggplot
data_long <- data_sorted_pert %>%
  select(sgn_number_factor, perturbation, 
         bing_height_log, lm_height_log, notlincs_height_log) %>%
  pivot_longer(
    cols = c(bing_height_log, lm_height_log, notlincs_height_log),
    names_to = "type",
    values_to = "height_log"
  ) %>%
  mutate(
    type = case_when(
      type == "bing_height_log" ~ "Best inferred",
      type == "lm_height_log" ~ "Landmark",
      type == "notlincs_height_log" ~ "Not LINCS"
    ),
    type = factor(type, levels = c("Not LINCS", "Best inferred", "Landmark"))
  )

# Цвета
colors <- c("Best inferred" = "steelblue",
            "Landmark" = "gray80", 
            "Not LINCS" = "black")

# Создаем график с geom_col
bp_sgnlen <- ggplot(data_long, 
                    aes(x = sgn_number_factor, 
                        y = height_log,
                        fill = type,
                        group = type)
) +
  geom_col(
    position = "stack",
    width = 0.74,
    color = "gray15",
    alpha = 0.9
  ) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(
    name = "Signatures"
  ) +
  scale_y_continuous(
    labels = function(x) {
      # Обратное преобразование для меток
      ifelse(x == 0, 0, round(10^(-x) - 1, 0))
    },
    breaks = function(limits) {
      # Создаем breaks на основе обратного преобразования
      breaks_abs <- c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
      breaks_transformed <- ifelse(breaks_abs == 0, 0, -log10(breaks_abs + 1))
      return(breaks_transformed[breaks_transformed >= limits[1] & breaks_transformed <= limits[2]])
    },
    # expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               vjust = 0.5,
                               size = 6),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(x = "Signatures", 
       y = "Number of genes per signature (log10 scale)", 
       title = "Signatures' composition")

print(bp_sgnlen)


# РИСУЕМ ДВА БАРПЛОТА ДЛЯ ОБЪЕДИНЕНИЯ

if (plots_mode == "dark") {
  upper_text_color = "white"
} else {
  upper_text_color = "black"
}

bp_top <- bp_strnum <- ggplot(data = data_sorted_pert,
                              mapping = aes(sgn_number_ordered, 
                                            as.numeric(method_res),
                                            fill = perturbation)) +
  geom_col(color = "gray15",
           width = 0.74) +
  theme_minimal() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  ) +
  geom_text(
    aes(label = label_symbol), # Указываем, какой столбец использовать как текст метки
    hjust = -0.09,            # Вертикальное выравнивание: -0.5 немного приподнимает текст над столбцом
    color = upper_text_color,         # Цвет текста
    size = 2.5,                 # Размер текста
    angle = 90,
    fontface = "bold"
  ) +
  labs(y = paste0("Number of perturbagen\nin ", method, " result")) +
  theme(
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               vjust = 0.5,
                               size = 6),
    legend.position = "top",
    legend.title = element_blank())
bp_top



bp_bottom <- ggplot(data_long, 
                    aes(x = sgn_number_factor, 
                        y = height_log,
                        fill = type,
                        group = type)
) +
  geom_col(
    position = "stack",
    width = 0.74,
    color = "gray15"
  ) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(
    name = "Signatures"
  ) +
  scale_y_continuous(
    labels = function(x) {
      # Обратное преобразование для меток
      ifelse(x == 0, 0, round(10^(-x) - 1, 0))
    },
    breaks = function(limits) {
      # Создаем breaks на основе обратного преобразования
      breaks_abs <- c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
      breaks_transformed <- ifelse(breaks_abs == 0, 0, -log10(breaks_abs + 1))
      return(breaks_transformed[breaks_transformed >= limits[1] & breaks_transformed <= limits[2]])
    },
    expand = expansion(mult = c(0.1, 0))
  ) +
  theme_minimal() +
  labs(y = "Amounts of genes\nper signature (log1p scale)") +
  theme(
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               vjust = 0.5,
                               size = 6),
    legend.position = "bottom",
    legend.title = element_blank()
  ) 

bp_bottom

# Верхний график: убираем нижние подписи / отступы
bp_top_adj <- bp_top +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 5, 0, 5)  # сверху, справа, снизу, слева
  ) + 
  coord_cartesian(ylim = c(0, 0.9)) # NA оставит нижнюю границу автоматической
bp_top_adj
# Нижний график: убираем верхние подписи / отступы
bp_bottom_adj <- bp_bottom +
  theme(
    plot.margin = margin(0, 5, 0, 5)
  )



colors <- c("Best inferred" = "steelblue",
            "Landmark" = "gray80", 
            "Not LINCS" = "black")


# === 0. ЗАДАЕМ ОБЩИЕ УРОВНИ ДЛЯ ОСИ X (ГАРАНТИЯ СОВПАДЕНИЯ) ===
common_x_levels <- levels(data_sorted_pert$sgn_number_ordered)

# 1. Верхний график
bp_top_adj <- bp_top +
  scale_fill_viridis_d(option = "magma", begin = 0.25) +
  # Жесткая привязка: одинаковые уровни и одинаковый отступ слева/справа
  scale_x_discrete(limits = common_x_levels, expand = expansion(add = 0.6)) +
  # coord_cartesian(..., expand = FALSE) устарел в ggplot2 >= 3.4.0. 
  # Используем scale_y_continuous для фиксации границ и отступов:
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

# 2. Нижний график
bp_bottom_adj <- bp_bottom +
  scale_fill_manual(values = colors) +
  # ИДЕНТИЧНЫЙ scale_x_discrete для пиксельного совпадения столбцов
  scale_x_discrete(limits = common_x_levels, expand = expansion(add = 0.6)) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0)), 
    labels = function(x) {
      ifelse(x == 0, 0, round(10^(-x) - 1, 0))
    },
    breaks = function(limits) {
      breaks_abs <- c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
      breaks_transformed <- ifelse(breaks_abs == 0, 0, -log10(breaks_abs + 1))
      return(breaks_transformed[breaks_transformed >= limits[1] & breaks_transformed <= limits[2]])
    }
  ) +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  )


final_barplot = ""

if (plots_mode == "dark") {
    # 3. Сборка с общими темными стилями и вертикальной легендой
    combined_plot_dark <- (bp_top_adj / bp_bottom_adj) + 
      plot_layout(heights = c(1, 1), guides = "collect") & 
      theme(
        # Темная тема
        panel.background = element_rect(fill = "grey30"),
        plot.background = element_rect(fill = "grey20", color = NA),
        panel.grid.major = element_line(color = "grey40"),
        panel.grid.minor = element_line(color = "grey35"),
        text = element_text(color = "white"),
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        axis.line = element_line(color = "grey60"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(color = "white"),
        
        # Склеивание панелей
        panel.spacing = unit(0, "lines"),
        
        # Вертикальная легенда
        legend.direction = "vertical",
        legend.box = "vertical"
      )
    
    # Визуализация
    print(combined_plot_dark)
    
    
    final_barplot = combined_plot_dark
    
    
    
    setwd(paste0(project_folder_name, inh_grp, "_inhibitors/", method, "/visualization"))
    
    ggsave(
      filename = paste0(method, "_dark_barplot_with_lengths.png"), # Имя файла и формат
      plot = combined_plot_dark,                                   # Какой объект графика сохранить
      width = 8,                                        # Ширина (по умолчанию в дюймах)
      height = 5,                                       # Высота (по умолчанию в дюймах)
      units = "in",                                     # Единицы измерения (in, cm, mm, px)
      dpi = 1000                                         # Разрешение (только для растровых форматов типа PNG, JPEG)
    )
    
    
} else {
  # РИСУЕМ ОБЪЕДИНЕННЫЙ БАПЛОТ
  combined_plot <- bp_top_adj / bp_bottom_adj + 
    plot_layout(heights = c(1, 1)) +
    theme(
      axis.text.x = element_text(angle = 90, 
                                 hjust = 1, 
                                 vjust = 0.5,
                                 size = 6)) +
    labs(x = "Signatures") +
    plot_annotation(
      theme = theme(plot.title = element_text(hjust = 0.5))) &
    theme(plot.margin = margin(0, 0, 0, 0)) 
  
  
  final_barplot = combined_plot
  
  
  setwd(paste0(project_folder_name, inh_grp, "_inhibitors/", method, "/visualization"))
  
  ggsave(
    filename = paste0(method, "_barplot_with_lengths.png"), # Имя файла и формат
    plot = combined_plot,                                   # Какой объект графика сохранить
    width = 8,                                        # Ширина (по умолчанию в дюймах)
    height = 5,                                       # Высота (по умолчанию в дюймах)
    units = "in",                                     # Единицы измерения (in, cm, mm, px)
    dpi = 1000                                         # Разрешение (только для растровых форматов типа PNG, JPEG)
  )
  
  
}




final_barplot







#str(data_sorted_pert$label_symbol)
#print(any(is.na(data_sorted_pert$label_symbol)))
#data_sorted_pert$label_symbol <- as.character(data_sorted_pert$label_symbol)






'''
# 1. Верхний график (оставляем magma или другой стиль)
bp_top_adj <- bp_top +
  scale_fill_viridis_d(option = "magma") + # Если здесь тоже нужны свои цвета, замените на scale_fill_manual
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  ) + 
  coord_cartesian(ylim = c(0, 1.1), expand = FALSE)

# 2. Нижний график (применяем вашу палитру напрямую)
bp_bottom_adj <- bp_bottom +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    # Добавляем небольшой отступ снизу (например, 5%)
    expand = expansion(mult = c(0.05, 0)), 
    labels = function(x) {
      ifelse(x == 0, 0, round(10^(-x) - 1, 0))
    },
    breaks = function(limits) {
      breaks_abs <- c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
      breaks_transformed <- ifelse(breaks_abs == 0, 0, -log10(breaks_abs + 1))
      return(breaks_transformed[breaks_transformed >= limits[1] & breaks_transformed <= limits[2]])
    }
  ) +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  )

# 3. Сборка с общими темными стилями и вертикальной легендой
combined_plot_dark <- (bp_top_adj / bp_bottom_adj) + 
  plot_layout(heights = c(1, 1), guides = "collect") & 
  theme(
    # Темная тема
    panel.background = element_rect(fill = "grey30"),
    plot.background = element_rect(fill = "grey20", color = NA),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_line(color = "grey35"),
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.line = element_line(color = "grey60"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color = "white"),
    
    # Склеивание
    panel.spacing = unit(0, "lines"),
    
    # Вертикальная легенда
    legend.direction = "vertical",
    legend.box = "vertical"
  )
'''



'''
# 1. Настраиваем ВЕРХНИЙ график
bp_top_adj <- bp_top +
  scale_fill_viridis_d(option = "magma") +
  # Важно: используем ту же переменную для X, что и в нижнем
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  scale_y_continuous(
    limits = c(0, 1.1),           # Четкая граница на 1.1
    expand = c(0, 0)
  ) +  
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 0, l = 5)
  )

# 2. Настраиваем НИЖНИЙ график
bp_bottom_adj <- bp_bottom +
  scale_fill_manual(values = colors) +
  # Важно: идентичный expand, чтобы столбцы не гуляли
  scale_x_discrete(name = "Signatures", expand = expansion(add = 0.6)) +
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0)), # Отступ снизу, чтобы не липло к оси
    labels = function(x) ifelse(x == 0, 0, round(10^(-x) - 1, 0)),
    breaks = function(limits) {
      breaks_abs <- c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
      breaks_transformed <- ifelse(breaks_abs == 0, 0, -log10(breaks_abs + 1))
      return(breaks_transformed[breaks_transformed >= limits[1] & breaks_transformed <= limits[2]])
    }
  ) +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
  )

# 3. Сборка через patchwork
combined_plot_dark <- (bp_top_adj / bp_bottom_adj) + 
  plot_layout(
    heights = c(1, 1), 
    guides = "collect"
  ) & 
  theme(
    # Темная тема
    panel.background = element_rect(fill = "grey30"),
    plot.background = element_rect(fill = "grey20", color = NA),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_line(color = "grey35"),
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.line = element_line(color = "grey60"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color = "white"),
    
    # ИДЕАЛЬНОЕ СКЛЕИВАНИЕ
    panel.spacing = unit(0, "lines"),
    legend.direction = "vertical",
    legend.box = "vertical"
  )

# Принудительное выравнивание панелей (align)
combined_plot_dark <- combined_plot + plot_layout(guides = 'collect')
'''


