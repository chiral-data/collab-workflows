#!/bin/bash
#

REGISTRY_URL="chiral.sakuracr.jp"
FORCE_REBUILD=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force|-f)
            FORCE_REBUILD=true
            shift
            ;;
        *)
            DIRECTORY_ARG="$1"
            shift
            ;;
    esac
done

function build_and_push_docker_image() {
    if [ -z "$1" ]; then
        echo "Error: No directory provided. Usage: build_and_push_docker_image <directory_path>"
        return 1
    fi

    local DIR="$1"

    echo "Entering $DIR"
    cd "$DIR" || {
        echo "Failed to enter $DIR"
        continue
    }

    fullname=$(basename "$PWD")
    app_name=$(echo "$fullname" | cut -d'_' -f1)
    date=$(echo "$fullname" | cut -d'_' -f2-)

    echo "Build container image for app $app_name, version: $date"
    image_name=${app_name}_${date}
    image_name_remote=${REGISTRY_URL}/${app_name}:${date}

    if [ "$FORCE_REBUILD" = true ]; then
        echo "Force rebuilding container image $image_name"
    elif docker manifest inspect "$image_name_remote" >/dev/null 2>&1; then
        echo "Container image $image_name_remote exist ... Skip build"
        return 0
    else
        echo "Building container image $image_name"
    fi

    echo ""
    docker build -t ${image_name} --platform linux/amd64 .
    docker run --gpus all ${image_name}

    echo "Tag container image $image_name"
    docker tag ${image_name} ${image_name_remote}

    echo "Push container image $image_name"
    docker push ${image_name_remote}

    cd - >/dev/null
}

if [ -n "$DIRECTORY_ARG" ]; then
    echo "Build the image for folder: $DIRECTORY_ARG"
    build_and_push_docker_image "$DIRECTORY_ARG"
else
    echo "Build all the images"

    SUBFOLDERS=$(find . -mindepth 2 -maxdepth 2 -type d ! -path "*/.git*")

    for DIR in $SUBFOLDERS; do
        build_and_push_docker_image "$DIR"
    done
fi
