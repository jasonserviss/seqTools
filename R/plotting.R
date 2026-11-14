#' coloursFromTargets
#'
#' Diverging color palette with 40 colors.
#'
#' @name coloursFromTargets
#' @rdname coloursFromTargets
#' @author Jason T. Serviss
#' @param pal character; A colour palette with length = length(markers).
#' @param counts matrix; A matrix containing counts.
#' @param markers character; The markers to evaluate. Must be present in
#'  rownames(counts).
#' @param ... additional arguments to pass on.
#' @keywords coloursFromTargets
NULL

#' @rdname coloursFromTargets
#' @importFrom dplyr group_by ungroup mutate arrange summarize select if_else n
#' @importFrom rlang .data
#' @importFrom tibble tibble add_column
#' @importFrom tidyr gather unnest spread
#' @importFrom purrr pmap pmap_chr
#' @importFrom grDevices col2rgb rgb
#' @importFrom readr parse_factor

coloursFromTargets <- function(
  pal, counts, markers, ...
){
  
  if(is.null(markers) | is.null(pal) | length(markers) == 1) {
    return(tibble('Sample' = colnames(counts)))
  }
  rs <- rowSums(counts[markers, ])
  if(any(rs == 0)) {
    g <- paste(markers[which(rs == 0)], collapse = ", ")
    mess <- paste0("The following genes have 0 counts for all samples ", g)
    stop(mess)
  }
  
  markers <- sort(markers)
  pal <- pal[1:length(markers)]
  
  counts[rownames(counts) %in% markers, ] |>
  matrix_to_tibble(data = _, rowname = 'geneName', drop = FALSE) |>
  pivot_longer(
    cols = !.data$geneName, names_to = 'Sample', values_to = 'count'
  ) |>
  #normalize
  group_by(.data$geneName) |>
  mutate('normalized' = normalize.vector(.data$count)) |>
  ungroup() |>
  #calculate fraction
  group_by(.data$Sample) |>
  mutate('fraction' = .data$normalized / sum(.data$normalized)) |>
  mutate(
    'fraction' = if_else(is.nan(.data$fraction), 1 / n(), .data$fraction)
  ) |>
  #setup initial hex colors
  arrange(.data$geneName) |>
  mutate('colint' = pal) |>
  ungroup() |>
  #convert to rgb and calculate new colors
  mutate('rgb' = pmap(
    list(.data$colint, .data$normalized, .data$fraction),
    function(x, y, z)  {
     mat <- as.vector((255 - ((255 - col2rgb(x)) * y)) * z)
     names(mat) <- c("r", "g", "b")
     mat
    }
  )) |>
  unnest_longer(col = c(.data$rgb), values_to = "rgb", indices_to = "col") |>
  group_by(.data$Sample, .data$col) |>
  summarize('sumRGB' = sum(.data$rgb) / 256, .groups = "drop_last") |>
  ungroup() |>
  pivot_wider(names_from=col, values_from = sumRGB) |>
  #convert back to hex
  mutate('Colour' = pmap_chr(
    list(.data$r, .data$g, .data$b),
    function(x, y, z) rgb(red = x, green = y, blue = z)
  )) |>
  select(-(.data$b:.data$r)) |>
  #fix factor levels so ggplot legend will cooperate
  #https://community.rstudio.com/t/drop-false-with-scale-fill-identity/5163/2
  mutate('Colour' = parse_factor(
    .data$Colour,
    levels = unique(c(.data$Colour, pal[!pal %in% .data$Colour]))
  ))
}

#' col40
#'
#' Diverging color palette with 40 colors.
#'
#' @name col40
#' @rdname col40
#' @author Jason T. Serviss
#' @keywords col40
#' @examples
#'
#' col40()
#'
#' @export
NULL

col40 <- function() {
  c(
  "#1c54a8", "#f63b32", "#00C2A0", "#FFDBE5", "#e0c48c", "#63FFAC", "#663000",
  "#e0a81c", "#385438", "#609060", "#6A3D9A", "#548495", "#A30059", "#8FB0FF",
  "#997D87", "#4FC601", "#8ca8c4", "#3B5DFF", "#BA0900", "#DDEFFF", "#7B4F4B",
  "#A1C299", "#0AA6D8", "#00846F", "#FFB500", "#C2FFED", "#A079BF", "#C0B9B2",
  "#C2FF99", "#00489C", "#6F0062", "#EEC3FF", "#922329", "#FFF69F", "#FF8A9A",
  "#B05B6F", "#7900D7", "#BC23FF", "#9B9700", "#0089A3"
  )
}
