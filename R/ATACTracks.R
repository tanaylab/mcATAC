#' ATACTracks objects
#'
#'
#' @description ATACTracks is a shallow object holding ATAC data over cells/metacells. Minimally it should include a vector of track names which hold the names of the tracks with the ATAC data.
#' McTracks extend the ATACTracks object by adding metadata and additional slots.
#'
#' @slot tracks vector of track names with the ATAC data.
#' @slot total_cov vector with total coverage of the tracks.
#' @slot marginal_track name of a track with marginal coverage.
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations. The constructor can also include or the name of a delimited file which contains such annotations.
#' @slot resolution the resolution of the tracks.
#' @slot window_size The size of the window used to smooth the counts in each track. If no smoothing was done - this should be equal to resolution/2
#'
#' @exportClass ATACTracks
ATACTracks <- setClass(
    "ATACTracks",
    slots = c(
        tracks = "character",
        total_cov = "numeric",
        marginal_track = "character_or_null",
        resolution = "numeric",
        window_size = "numeric"
    ),
    contains = c("ATAC", "VIRTUAL")
)

setMethod(
    "initialize",
    signature = "ATACTracks",
    definition = function(.Object, tracks, genome, id = NULL, description = NULL, path = "", marginal_track = NULL, resolution = NULL, window_size = resolution / 2, order = NULL) {
        .Object <- make_atac_tracks_object(.Object, tracks, genome, id, description, path = path, marginal_track = marginal_track, resolution = resolution, window_size = window_size, order = order)
        return(.Object)
    }
)

make_atac_tracks_object <- function(obj, tracks, genome, id, description, path = "", total_cov = NULL, marginal_track = NULL, resolution = NULL, window_size = resolution / 2, order = NULL) {
    obj <- make_atac_object(obj, genome, id, description, path)

    walk(tracks, ~ {
        if (!gtrack.exists(.x)) {
            cli_abort("Track {.x} does not exist")
        }
    })

    if (is.null(total_cov)) {
        total_cov <- plyr::laply(tracks, get_track_total_coverage, .parallel = getOption("mcatac.parallel"))
    }

    if (is.null(resolution)) {
        resolution <- gtrack.info(tracks[1])$bin.size
    }

    if (is.null(window_size)) {
        window_size <- resolution
    }

    if (is.null(order)) {
        order <- 1:length(tracks)
    } else {
        if (length(order) != length(tracks)) {
            cli_abort("Number of columns in the matrix is not equal to the number of tracks.")
        }
        if (any(order < 1) || any(order > length(tracks))) {
            cli_abort("Order vector contains values outside the range of the tracks.")
        }
    }

    obj@tracks <- tracks
    obj@total_cov <- total_cov
    obj@marginal_track <- marginal_track
    obj@resolution <- resolution
    obj@window_size <- window_size
    obj@order <- order
    return(obj)
}

#' McTracks
#'
#' @rdname ATACTracks
#' @exportClass McTracks
McTracks <- setClass(
    "McTracks",
    slots = c(
        metacells = "character"
    ),
    contains = c("ATACTracks")
)

setMethod(
    "initialize",
    signature = "McTracks",
    definition = function(.Object, tracks, genome, metacells = NULL, id = NULL, description = NULL, path = "", marginal_track = NULL, resolution = NULL, window_size = resolution / 2, order = NULL) {
        if (is.null(metacells)) {
            mc_nums <- stringr::str_extract_all(string = tracks, pattern = "mc\\d*$") %>%
                unlist() %>%
                gsub("mc", "", .) %>%
                as.numeric()
            tracks <- tracks[order(mc_nums)]
            metacells <- mc_nums[order(mc_nums)]
        } else {
            if (length(metacells) != length(tracks)) {
                cli_abort("Number of tracks and metacells must be equal")
            }
        }
        .Object <- make_atac_tracks_object(.Object, tracks, genome, id, description, path = path, marginal_track = marginal_track, resolution = resolution, window_size = window_size, order = order)
        .Object@metacells <- as.character(metacells)

        return(.Object)
    }
)

#' Create an McTracks object
#'
#' @description Create an McTracks object either from a list of tracks (\code{tracks} parameter) or from a track prefix (\code{track_prefix} parameter). The tracks are assumed to be in the form of ".+_mc{metacell_number}", if this is not the case, the \code{metacells} parameter can be used to specify the metacell number for each track.
#'
#' @param genome genome assembly of the tracks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param tracks a vector of tracks to include in the object.
#' @param track_prefix prefix for tracks (instead of explicit track names at the \code{tracks} parameter). Ignored if \code{tracks} is not NULL. If
#' a track named "{track_prefix}.marginal" exists, it will be used as the marginal track.
#' @param metacells a vector of metacells to include in the object (optional). If not provided, the metacells are inferred from the track names using the pattern "{track_prefix}.mc\\d*$". In such a case, the tracks are ordered according to the metacell number.
#' @param marginal_track a track to use as the marginal track (optional). If not provided, and \code{track_prefix} was given, the track named "{track_prefix}.marginal" will be used if it exists.
#' @param id an identifier for the object, e.g. "pbmc".
#' @param path path to a directory containing the raw data. Usually inherited from the counts object.
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#' @param resolution the resolution of the tracks.
#' @param window_size The size of the window used to smooth the counts in each track. If no smoothing was done - this should be equal to resolution / 2.
#' @param order the order of the tracks in the object. If not provided, the tracks are ordered according to the metacell number.
#'
#' @return a McTracks object.
#'
#' @examples
#' \dontrun{
#' mct_create(track_prefix = "pbmc_mc", genome = "hg38", id = "pbmc", description = "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)")
#' }
#'
#' @export
mct_create <- function(genome, tracks = NULL, track_prefix = NULL, metacells = NULL, marginal_track = NULL, id = NULL, description = NULL, path = NULL, metadata = NULL, resolution = NULL, window_size = resolution / 2, order = NULL) {
    gset_genome(genome)
    if (is.null(tracks)) {
        if (is.null(track_prefix)) {
            cli_abort("Either {.field tracks} or {.field track_prefix} must be provided")
        }
        tracks <- gtrack.ls(track_prefix)
        marginal_track_pref <- paste0(track_prefix, ".marginal")
        tracks <- tracks[tracks != marginal_track_pref]
        marginal_track <- marginal_track %||% marginal_track_pref
    }

    if (!gtrack.exists(marginal_track)) {
        marginal_track <- NULL
    }

    path <- path %||% ""

    res <- new("McTracks", tracks = tracks, genome = genome, metacells = metacells, id = id, description = description, path = path, marginal_track = marginal_track, resolution = resolution, window_size = window_size, order = order)

    cli_alert_success("Created a new McTracks object with {.val {length(res@metacells)}} metacells.")
    return(res)
}

#' @export
#' @noRd
setMethod(
    "show",
    signature = "McTracks",
    definition = function(object) {
        print_atac_tracks_object(object, "McTracks", "metacell", "metacell")
    }
)

print_atac_tracks_object <- function(object, object_type, column_type, md_column) {
    cli::cli_text("{.cls {object_type}} object with {.val {length(object@metacells)}} {column_type}s from {.field {object@genome}}.")
    if (object@id != "") {
        cli::cli_text(c("id: {.val {object@id}}"))
    }
    if (object@description != "") {
        cli::cli_text(c("description: {.val {object@description}}"))
    }
    if (object@path != "") {
        cli::cli_text(c("Loaded from: {.file {object@path}}"))
    }
    cli::cli_text(c("resolution: {.val {object@resolution}}"))
    cli::cli_text(c("window_size: {.val {object@window_size}}"))

    cli::cli_text("Slots include:")
    cli_ul(c("{.code @tracks}: a vector with track names"))
    cli_ul(c("{.code @total_cov}: a vector with total coverage per track"))
    cli_ul(c("{.code @genome}: genome assembly of the peaks"))
    if (!is.null(object@marginal_track)) {
        cli_ul(c("{.code @marginal_track}: a track with marginal counts (sum of all the tracks)"))
    }
    if (object_type == "McTracks") {
        cli_ul(c("{.code @metacells}: a character vector with metacell ids"))
        if (has_rna(object)) {
            cli_ul(c("{.code @rna_egc}: a numeric matrix which contains normalized RNA expression per gene (rows) per metacell (columns)."))
        }
    }
    if (!is.null(object@metadata)) {
        cli_ul(c("{.code @metadata}: a tibble with a column called '{md_column}' and additional {column_type} annotations."))
    }
    cli::cli_text("Tracks (only first are shown):")
    walk(object@tracks[1:min(5, length(object@tracks))], function(x) {
        cli_ul(c("{.val {x}}"))
    })
}
