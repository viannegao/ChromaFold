import os
import numpy as np

import torch
import torch.nn as nn

from tqdm import tqdm
from datetime import datetime


def train_bulk_only(
    train_loader, model, criterion, optimizer, device, min_stripe_signal
):
    model.train()
    running_loss = 0

    with tqdm(train_loader, unit="batch") as tepoch:
        for X1, X1_rev, X2, X2_rev, y_true, y_true_rev in tepoch:
            optimizer.zero_grad()

            # left interactions
            X1 = X1.to(device)
            X2 = X2.to(device)
            X = torch.cat((X1, X2), axis=1)
            y_true = y_true.to(device)
            y_hat = model(X)
            loss = criterion(y_hat, y_true, min_stripe_signal)

            # right interactions
            X1_rev = X1_rev.to(device)
            X2_rev = X2_rev.to(device)
            X_rev = torch.cat((X1_rev, X2_rev), axis=1)
            y_true_rev = y_true_rev.to(device)
            y_hat_rev = model(X_rev)
            loss += criterion(y_hat_rev, y_true_rev, min_stripe_signal)

            running_loss += loss.item()

            # Backward pass
            loss.backward()
            optimizer.step()

    epoch_loss = running_loss / len(train_loader.dataset)

    return model, optimizer, epoch_loss


def validate_bulk_only(val_loader, model, criterion, device, min_stripe_signal):
    model.eval()
    running_loss = 0

    for X1, X1_rev, X2, X2_rev, y_true, y_true_rev in val_loader:

        # left interactions
        X1 = X1.to(device)
        X2 = X2.to(device)
        X = torch.cat((X1, X2), axis=1)
        y_true = y_true.to(device)
        y_hat = model(X)
        loss = criterion(y_hat, y_true, min_stripe_signal)

        # right interactions
        X1_rev = X1_rev.to(device)
        X2_rev = X2_rev.to(device)
        X_rev = torch.cat((X1_rev, X2_rev), axis=1)
        y_true_rev = y_true_rev.to(device)
        y_hat_rev = model(X_rev)
        loss += criterion(y_hat_rev, y_true_rev, min_stripe_signal)

        running_loss += loss.item()

    epoch_loss = running_loss / len(val_loader.dataset)

    return model, epoch_loss


def train(train_loader, model, criterion, optimizer, device, min_stripe_signal):
    model.train()
    running_loss = 0
    with tqdm(train_loader, unit="batch") as tepoch:
        for Z2, Z2_rev, X1, X1_rev, X2, X2_rev, y_true, y_true_rev in tepoch:
            optimizer.zero_grad()

            # left interactions
            X1 = X1.to(device)
            X2 = X2.to(device)
            X = torch.cat((X1, X2), axis=1)
            Z2 = Z2.to(device)
            y_true = y_true.to(device)

            y_hat = model(Z2, X)
            loss = criterion(y_hat, y_true, min_stripe_signal)

            # right interactions
            X1_rev = X1_rev.to(device)
            X2_rev = X2_rev.to(device)
            X_rev = torch.cat((X1_rev, X2_rev), axis=1)
            Z2_rev = Z2_rev.to(device)
            y_true_rev = y_true_rev.to(device)

            y_hat_rev = model(Z2_rev, X_rev)
            loss += criterion(y_hat_rev, y_true_rev, min_stripe_signal)

            running_loss += loss.item()  # * X.size(0)

            # Backward pass
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 100000)
            optimizer.step()

    epoch_loss = running_loss / len(train_loader.dataset)
    return model, optimizer, epoch_loss


def validate(val_loader, model, criterion, device, min_stripe_signal):
    """
    Function for the validation step of the training loop
    """
    model.eval()
    running_loss = 0

    for Z2, Z2_rev, X1, X1_rev, X2, X2_rev, y_true, y_true_rev in val_loader:
        # left interactions
        X1 = X1.to(device)
        X2 = X2.to(device)
        X = torch.cat((X1, X2), axis=1)
        Z2 = Z2.to(device)
        y_true = y_true.to(device)

        y_hat = model(Z2, X)
        loss = criterion(y_hat, y_true, min_stripe_signal)

        # right interactions
        X1_rev = X1_rev.to(device)
        X2_rev = X2_rev.to(device)
        X_rev = torch.cat((X1_rev, X2_rev), axis=1)
        Z2_rev = Z2_rev.to(device)
        y_true_rev = y_true_rev.to(device)

        y_hat_rev = model(Z2_rev, X_rev)
        loss += criterion(y_hat_rev, y_true_rev, min_stripe_signal)

        running_loss += loss.item()  # * X.size(0)

    epoch_loss = running_loss / len(val_loader.dataset)

    return model, epoch_loss


def training_loop(
    train,
    validate,
    model,
    ct_list,
    seed,
    criterion,
    optimizer,
    scheduler,
    train_loader,
    valid_loader,
    epochs,
    patience,
    save_path,
    device,
    min_stripe_signal,
    print_every=1,
):
    """
    Function defining the entire training loop
    """
    # set objects for storing metrics
    counter = 0
    best_loss = 1e10
    train_losses = []
    valid_losses = []

    SAVEPATH = os.path.join(
        save_path, "{}_CTCFmotifScore_seed{}.pth.tar".format("_".join(ct_list), seed)
    )

    # Train model
    for epoch in range(epochs):

        # training
        model, optimizer, train_loss = train(
            train_loader, model, criterion, optimizer, device, min_stripe_signal
        )
        train_losses.append(train_loss)

        # validation
        with torch.no_grad():
            model, valid_loss = validate(
                valid_loader, model, criterion, device, min_stripe_signal
            )
            valid_losses.append(valid_loss)

        counter += 1

        if valid_loss < best_loss:
            counter = 0
            best_loss = valid_loss
            print("++++saving+++++")
            torch.save(model.state_dict(), SAVEPATH)

        if epoch % print_every == (print_every - 1):
            print(
                f"{datetime.now().time().replace(microsecond=0)} --- "
                f"Epoch: {epoch}\t"
                f"Train loss: {train_loss:.4f}\t"
                f"Valid loss: {valid_loss:.4f}\t"
            )

        if counter > patience:
            print("Early stopping.")
            break

        scheduler.step()

    return model, optimizer, (train_losses, valid_losses)


def weighted_mse_loss(input, target, min_stripe_signal=-99999):
    # do not compute loss on low-mappability regions
    weight = target.sum(1) > min_stripe_signal
    loss = (input - target) ** 2
    loss = (weight * loss.T).T
    loss = torch.sum(loss)

    return loss
